#!/bin/bash
# mutect 0.0.1

main() {
    set -ex -o pipefail
    # which bcftools
    # if [[ -f /usr/bin/bcftools && -f /usr/bin/tabix ]]; then
    #     /usr/bin/bcftools -v
    #     /usr/bin/tabix --version
    # fi
    if [[ ! -f /usr/local/bin/Mutect2.sh ]]; then
        echo 'File "/usr/local/bin/Mutect2.sh" is not there, aborting.'
        exit
    fi
    # dx download UKBB_Exome_2021:CH_Exome/Inputs/GRCh38_full_analysis_set_plus_decoy_hla.fa  -o /usr/local/reference/
    # dx download UKBB_Exome_2021:CH_Exome/Inputs/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai  -o /usr/local/reference/
    # dx download UKBB_Exome_2021:CH_Exome/Inputs/GRCh38_full_analysis_set_plus_decoy_hla.dict  -o /usr/local/reference/
    # if [[ ! -f UKBB_Exome_2021:CH_Exome/Inputs/GRCh38_full_analysis_set_plus_decoy_hla.fa ]]; then
    #     echo 'File "UKBB_Exome_2021:/CH_Exome/Inputs/GRCh38_full_analysis_set_plus_decoy_hla.fa" is not there, aborting.'
    #     exit
    # fi
    
    # ls 
    # ls /
    # if [[ ! -f /gatk/gatk-4.2.1.0/gatk ]]; then
    #     echo 'File "/gatk/gatk-4.2.1.0/gatk" is not there, aborting.'
    #     exit
    # else
    #     /gatk/gatk-4.2.1.0/gatk Tools -h
    # fi
    #dx ls /CH_Exome/Inputs
    #dx ls UKBB_Exome_2021:/CH_Exome/Inputs/GRCh38_full_analysis_set_plus_decoy_hla.fa*

    # echo "Value of reference: '$reference'"
    # echo "Value of reference: '$reference_index'"
    # echo "Value of reference: '$reference_dict'"
    # echo "Value of tumor_bam: '$tumor_bam'"
    # echo "Value of tumor_bam: '$tumor_bai'"
    # echo "Value of normal_bam: '$normal_bam'"
    # echo "Value of normal_bam: '$normal_bai'"
    echo "Value of interval_list: '$interval_list'"
    # echo "Value of interval_lists_22: '$interval_lists_22'"
    # echo "Value of interval_lists_22[@]: ${interval_lists_22[@]}"
    #echo "Value of interval_lists_22: '${interval_lists_22}'"
    echo "Value of dockerimage_picard: '$dockerimage_picard'"
    # echo "Value of dockerimage_picard_path: '$dockerimage_picard_path'"
    # echo "Value of dockerimage_gatk: '$dockerimage_gatk'"
    # echo "Value of dockerimage_gatk_path: '$dockerimage_gatk_path'"
    # echo "Value of dockerimage_bcftools: '$dockerimage_bcftools'"
    # echo "Value of dockerimage_bcftools_path: '$dockerimage_bcftools_path'"
    
    #echo "${interval_lists_22[@]/#/-i array_of_scattered_input:array:file=}"
    #echo "array_of_scattered_input:array:file=${interval_lists_22[@]}"
    # scatter_job=$(dx-jobutil-new-job scatter \
    #     -i "input_to_scatter=${interval_list}" \
    #     -i "dockerimage_picard=${dockerimage_picard}")
    # scatter_job=$(dx-jobutil-new-job scatter \
    #     -i "input_to_scatter=${interval_list}")
    # dx wait ${scatter_job}
    # test=""
    # for file_obj in ${interval_lists_22[@]}
    # do
    #     test="$test -i array_of_scattered_input:array:file=${file_obj}"
    # done
    #test2="${interval_lists_22[@]/#/}"
    interval_files=()
    for file in "${interval_lists_22[@]}"; do
        interval_files+=("-iarray_of_scattered_input=${file}")
    done
    # echo "Value of files: '${files[@]}'"
    #-i array_of_scattered_input:array:jobref=${scatter_job}:array_of_scattered_input \
    #-i array_of_scattered_input:array:file="${interval_lists_22[@]}" \
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

    
    output_name='merged.mutect.filtered.vcf.gz'
    postprocess_job=$(dx-jobutil-new-job postprocess \
        -i process_outputs:array:jobref=${map_job}:process_outputs \
        -i output_name:string=${output_name} \
        --depends-on $map_job)
        # -i "dockerimage_bcftools=${dockerimage_bcftools}" \

    dx-jobutil-add-output vcf --class=jobref "$postprocess_job":vcf 
    dx-jobutil-add-output vcf_index --class=jobref "$postprocess_job":vcf_index
}

scatter() {
    set -ex -o pipefail
    echo "Value of input_to_scatter: '${input_to_scatter}'"
    # echo "Value of dockerimage_picard: '${dockerimage_picard}'"
    # echo "Value of dockerimage_picard_path: '${dockerimage_picard_path}'"
    echo "scattering"
    dx-download-all-inputs --parallel
  
    # Fill in code here to do whatever is necessary to scatter the
    # input.
    #docker load -i /picard.tar.gz
    #docker load -i ${dockerimage_picard_path}
    #which java
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
    # echo "Value of reference: '${reference}'"
    # echo "Value of reference_index: '${reference_index}'"
    # echo "Value of tumor_bai: '${tumor_bai}'"
    # echo "Value of normal_bai: '${normal_bai}'"
    # echo "Value of dockerimage_gatk: '$dockerimage_gatk'"
    # echo "Value of dockerimage_gatk_path: '$dockerimage_gatk_path'"

    ## might need counter i.e. a "shard" like Cromwell
    #dx download "$dockerimage_gatk"
    #docker load -i "$dockerimage_gatk_name"
    process_jobs=()
    i=1
    for interval_file in "${array_of_scattered_input[@]}"
    do
        echo "Value of interval_file: '${interval_file}'"
        
        process_jobs+=($(dx-jobutil-new-job process -i "interval_file=${interval_file}" -i "reference=${reference}" -i "reference_index=${reference_index}" -i "reference_dict=${reference_dict}" -i "tumor_bam=${tumor_bam}" -i "tumor_bai=${tumor_bai}" -i "normal_bam=${normal_bam}" -i "normal_bai=${normal_bai}" -i shard:string=shard-$i -i "dockerimage_gatk=${dockerimage_gatk}"))
        # process_jobs+=($(dx-jobutil-new-job process -i "interval_file=${interval_file}" -i "reference=${reference}" -i "reference_index=${reference_index}" -i "reference_dict=${reference_dict}" -i "tumor_bam=${tumor_bam}" -i "tumor_bai=${tumor_bai}" -i "normal_bam=${normal_bam}" -i "normal_bai=${normal_bai}" -i shard:string=shard-$i))
        echo "${process_jobs[@]}"
        i=$((i+1))
        echo "$i.after"
    done

    for process_job in "${process_jobs[@]}"
    do
        dx-jobutil-add-output process_outputs --class=array:jobref "${process_job}:vcf" --array
        # dx-jobutil-add-output process_outputs --class=file $(dx upload --brief "${process_job}:vcf") --array
        dx-jobutil-add-output process_outputs --class=array:jobref "${process_job}:vcf_index" --array 
        # dx-jobutil-add-output process_outputs --class=file $(dx upload --brief "${process_job}:vcf_index") --array
    done

}

process() {
    set -ex -o pipefail
   
    # echo "processing"
    # echo "Value of interval_file: '${interval_file}'"
    # echo "Value of reference: '${reference}'"
    # echo "Value of reference_index: '${reference_index}'"
    # echo "Value of reference_dict: '${reference_dict}'"
    # echo "Value of tumor_bam: '${tumor_bam}'"
    # echo "Value of tumor_bai: '${tumor_bai}'"
    # echo "Value of normal_bam: '${normal_bam}'"
    # echo "Value of normal_bai: '${normal_bai}'"
    # echo "Value of shard: '${shard}'"
    # echo "Value of dockerimage_gatk: '$dockerimage_gatk'"
    # echo "Value of dockerimage_gatk_path: '$dockerimage_gatk_path'"
    
    # Fill in code here to process the input and create output.
    dx-download-all-inputs --parallel
    mv $tumor_bai_path ~/in/tumor_bam
    mv $normal_bai_path ~/in/normal_bam
    mv $reference_index_path ~/in/reference
    mv $reference_dict_path ~/in/reference
    #docker load -i /gatk.tar.gz
    docker load -i ${dockerimage_gatk_path}
  
    #${reference_path}
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin/:/usr/local/bin -w /home/dnanexus broadinstitute/gatk:4.2.1.0 \
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

    dx-jobutil-add-output vcf --class=file "${output_name}"  
    dx-jobutil-add-output vcf_index --class=file "${output_name}.tbi" 
}
