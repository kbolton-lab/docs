#!/bin/bash
# mutect 0.0.1

main() {
    set -ex -o pipefail
    dx download "$tumor_bam"
    if [[ ! -f /usr/local/bin/Mutect2.sh ]]; then
        echo 'File "/usr/local/bin/Mutect2.sh" is not there, aborting.'
        exit
    fi
    # if [[ ! -f /usr/local/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa ||
    #       ! -f /usr/local/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai ||
    #       ! -f /usr/local/reference/GRCh38_full_analysis_set_plus_decoy_hla.dict ]]; then
    #     echo 'File "/usr/local/bin/Mutect2.sh" is not there, aborting.'
    #     exit
    # fi
    echo "Value of output_name: '$output_name'"
    new_name=$(/usr/bin/samtools view -H $tumor_bam_name | /usr/bin/perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)
    echo $new_name
   
    
    echo "Value of interval_list: '$interval_list'"
    echo "Value of dockerimage_picard: '$dockerimage_picard'"
   
    interval_files=()
    for file in "${interval_lists_22[@]}"; do
        interval_files+=("-iarray_of_scattered_input=${file}")
    done

    # if using scatter()
    #-i array_of_scattered_input:array:jobref=${scatter_job}:array_of_scattered_input \
 
    map_job=$(dx-jobutil-new-job map \
        "${interval_files[@]}" \
        -i "reference=${reference}" \
        -i "reference_index=${reference_index}" \
        -i "reference_dict=${reference_dict}" \
        -i "tumor_bam=${tumor_bam}" \
        -i "tumor_bai=${tumor_bai}" \
        -i "normal_bam=${normal_bam}" \
        -i "normal_bai=${normal_bai}" \
        -i "dockerimage_gatk=${dockerimage_gatk}")
        # -i "dockerimage_gatk=${dockerimage_gatk}"

    
    # output_name='merged.mutect.filtered.vcf.gz'
    postprocess_job=$(dx-jobutil-new-job postprocess \
        -i process_outputs:array:jobref=${map_job}:process_outputs \
        -i output_name:string=${output_name} \
        --depends-on $map_job)
        # -i "dockerimage_bcftools=${dockerimage_bcftools}" \

    dx-jobutil-add-output vcf --class=jobref "$postprocess_job":vcf 
    dx-jobutil-add-output vcf_index --class=jobref "$postprocess_job":vcf_index
    dx-jobutil-add-output tumor_sample_name --class=string "$new_name"
}

scatter() {
    set -ex -o pipefail
    echo "Value of input_to_scatter: '${input_to_scatter}'"

    echo "scattering"
    dx-download-all-inputs --parallel
  
    output_dir="out"
    mkdir $output_dir
    scatter_count=10

    # docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/bin/:/usr/local/bin -w /home/dnanexus broadinstitute/picard:2.23.6 \
        # /usr/bin/perl /usr/local/bin/split_interval_list_helper.pl /home/dnanexus/${output_dir} ${input_to_scatter_path} ${scatter_count}
    /usr/bin/perl /usr/local/bin/split_interval_list_helper.pl /home/dnanexus/${output_dir} ${input_to_scatter_path} ${scatter_count}

    scattered_input=( $(ls ${output_dir}/*.interval_list) )
    
    for piece in "${scattered_input[@]}"
    do
        dx-jobutil-add-output array_of_scattered_input $(dx upload --brief $piece) --class=array:file --array
    done
   
}

map() {
    set -ex -o pipefail
    
    echo "Value of array_of_scattered_input: '${array_of_scattered_input[@]}'" 

    ## might need counter i.e. a "shard" like Cromwell
    #dx download "$dockerimage_gatk"
    #docker load -i "$dockerimage_gatk_name"
    process_jobs=()
    i=1
    for interval_file in "${array_of_scattered_input[@]}"
    do
        echo "Value of interval_file: '${interval_file}'"
        
        process_jobs+=($(dx-jobutil-new-job process -i "interval_file=${interval_file}" -i "reference=${reference}" -i "reference_index=${reference_index}" -i "reference_dict=${reference_dict}" -i "tumor_bam=${tumor_bam}" -i "tumor_bai=${tumor_bai}" -i "normal_bam=${normal_bam}" -i "normal_bai=${normal_bai}" -i shard:string=shard-$i -i "dockerimage_gatk=${dockerimage_gatk}"))
        # process_jobs+=($(dx-jobutil-new-job process -i "interval_file=${interval_file}" -i "tumor_bam=${tumor_bam}" -i "tumor_bai=${tumor_bai}" -i "normal_bam=${normal_bam}" -i "normal_bai=${normal_bai}" -i shard:string=shard-$i -i "dockerimage_gatk=${dockerimage_gatk}"))
        # process_jobs+=($(dx-jobutil-new-job process -i "interval_file=${interval_file}" -i "reference=${reference}" -i "reference_index=${reference_index}" -i "reference_dict=${reference_dict}" -i "tumor_bam=${tumor_bam}" -i "tumor_bai=${tumor_bai}" -i "normal_bam=${normal_bam}" -i "normal_bai=${normal_bai}" -i shard:string=shard-$i))
        echo "${process_jobs[@]}"
        i=$((i+1))
        echo "$i.after"
    done

    for process_job in "${process_jobs[@]}"
    do
        dx-jobutil-add-output process_outputs --class=array:jobref "${process_job}:vcf" --array
        dx-jobutil-add-output process_outputs --class=array:jobref "${process_job}:vcf_index" --array 
    done

}

process() {
    set -ex -o pipefail
    # if [[ ! -f /usr/local/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa ||
    #       ! -f /usr/local/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai ||
    #       ! -f /usr/local/reference/GRCh38_full_analysis_set_plus_decoy_hla.dict ]]; then
    #     echo 'File "/usr/local/bin/Mutect2.sh" is not there, aborting.'
    #     exit
    # fi
    # Fill in code here to process the input and create output.
    dx-download-all-inputs --parallel
    mv $tumor_bai_path ~/in/tumor_bam
    mv $normal_bai_path ~/in/normal_bam
    mv $reference_index_path ~/in/reference
    mv $reference_dict_path ~/in/reference
    #docker load -i /gatk.tar.gz
    docker load -i ${dockerimage_gatk_path}
  
    #${reference_path}
    #/usr/local/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin/:/usr/local/bin -v /usr/local/reference/:/usr/local/reference/ -w /home/dnanexus broadinstitute/gatk:4.2.1.0 \
        /bin/bash /usr/local/bin/Mutect2.sh ${shard}.mutect.vcf.gz ${reference_path} ${tumor_bam_path} ${normal_bam_path} ${interval_file_path} ${shard}.filtered.mutect.vcf.gz
    # /bin/bash /usr/local/bin/Mutect2.sh ${shard}.mutect.vcf.gz ${reference_path} ${tumor_bam_path} ${normal_bam_path} ${interval_file_path}
    

    vcf=$(dx upload ${shard}.filtered.mutect.vcf.gz --brief)
    dx-jobutil-add-output vcf --class=file "$vcf" 

    vcf_index=$(dx upload ${shard}.filtered.mutect.vcf.gz.tbi --brief)
    dx-jobutil-add-output vcf_index --class=file "$vcf_index"
}

postprocess() {
    set -ex -o pipefail
    if [[ -f /usr/bin/bcftools && -f /usr/bin/tabix ]]; then
        /usr/bin/bcftools -v
        /usr/bin/tabix --version
    else 
        exit
    fi
    echo "Value of process_outputs: '${process_outputs[@]}'"
    echo "Value of output_name: '${output_name}'"
    # echo "Value of dockerimage_bcftools: '$dockerimage_bcftools'"

    # download vcfs and vcf indices
    mkdir -p /tmp/mutect
    cd /tmp/mutect
    for process_output in "${process_outputs[@]}" 
    do
        dx download "$process_output"
        echo "$process_output_path"
    done
    ls -1sh /tmp/mutect/*.vcf.gz*

    cd /home/dnanexus
    #docker load -i /bcftools.tar.gz
    #dx download ${dockerimage_bcftools}
    #docker load -i ${dockerimage_bcftools_path}

    # docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/bin/:/usr/local/bin -w /home/dnanexus kboltonlab/sam_bcftools_tabix_bgzip:1.0 \
    #     /usr/local/bin/bcftools merge --merge none -Oz -o ${output_name} /tmp/mutect/*.vcf.gz && tabix ${output_name}
    /usr/bin/bcftools concat --allow-overlaps --remove-duplicates -Oz -o ${output_name} --threads 8 /tmp/mutect/*.vcf.gz && /usr/bin/tabix ${output_name}

    vcf=$(dx upload "${output_name}" --brief)
    dx-jobutil-add-output vcf --class=file "$vcf" 

    vcf_index=$(dx upload "${output_name}.tbi" --brief)
    dx-jobutil-add-output vcf_index --class=file "$vcf_index"

}
