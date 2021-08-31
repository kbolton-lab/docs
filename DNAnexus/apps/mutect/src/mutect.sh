#!/bin/bash
# mutect 0.0.1

main() {
    set -ex -o pipefail
    echo "8/24/2020"


    dx download "$tumor_bam"
    if [[ ! -f /usr/local/bin/Mutect2.sh ]]; then
        echo 'File "/usr/local/bin/Mutect2.sh" is not there, aborting.'
        exit
    fi
  
    #echo "Value of output_name: '$output_name'"
    ## for extracting tumor
    tumor_sample_name=$(/usr/bin/samtools view -H $tumor_bam_name | /usr/bin/perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)
    #echo $tumor_sample_name
   
    echo "Value of interval_list: '$interval_list'"
    echo "Value of dockerimage_picard: '$dockerimage_picard'"
   
    interval_files=()
    for file in "${interval_lists_array[@]}"; do
        interval_files+=("-iarray_of_scattered_input=${file}")
    done

 
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

    ## get <eid>_23153_0_0 from <eid>_23153_0_0.bqsr.bam
    eid_nameroot=$(echo $tumor_bam_name | cut -d'.' -f1)
    output_name="$eid_nameroot.normalized.merged.mutect.filtered.vcf.gz"
    ## add reference to turn on bcftools norm left-align
    postprocess_job=$(dx-jobutil-new-job postprocess \
        -i process_outputs:array:jobref=${map_job}:process_outputs \
        -i tumor_sample_name:string=${tumor_sample_name} \
        -i output_name:string=${output_name} \
        -i "reference:file=${reference}" \
        -i "reference_index:file=${reference_index}" \
        -i "tumor_bam=${tumor_bam}" \
        -i "tumor_bai=${tumor_bai}" \
        -i "dockerimage_gatk=${dockerimage_gatk}" \
        --depends-on $map_job)

    dx-jobutil-add-output vcf --class=jobref "$postprocess_job":vcf 
    dx-jobutil-add-output vcf_index --class=jobref "$postprocess_job":vcf_index
    dx-jobutil-add-output tumor_sample_name --class=string "$tumor_sample_name"
}

# scatter() {
#     set -ex -o pipefail
#     echo "Value of input_to_scatter: '${input_to_scatter}'"
#     dx-download-all-inputs --parallel
#     output_dir="out"
#     mkdir $output_dir
#     scatter_count=10

#     docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/bin/:/usr/local/bin -w /home/dnanexus broadinstitute/picard:2.23.6 \
#         /usr/bin/perl /usr/local/bin/split_interval_list_helper.pl /home/dnanexus/${output_dir} ${input_to_scatter_path} ${scatter_count}

#     scattered_input=( $(ls ${output_dir}/*.interval_list) )
    
#     for piece in "${scattered_input[@]}"
#     do
#         dx-jobutil-add-output array_of_scattered_input $(dx upload --brief $piece) --class=array:file --array
#     done
# }

map() {
    set -ex -o pipefail
    
    echo "Value of array_of_scattered_input: '${array_of_scattered_input[@]}'" 

    process_jobs=()
    i=1
    for interval_file in "${array_of_scattered_input[@]}"
    do
        echo "Value of interval_file: '${interval_file}'"
        
        process_jobs+=($(dx-jobutil-new-job process -i "interval_file=${interval_file}" -i "reference=${reference}" -i "reference_index=${reference_index}" -i "reference_dict=${reference_dict}" -i "tumor_bam=${tumor_bam}" -i "tumor_bai=${tumor_bai}" -i "normal_bam=${normal_bam}" -i "normal_bai=${normal_bai}" -i shard:string=shard-$i -i "dockerimage_gatk=${dockerimage_gatk}"))

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

    dx-download-all-inputs --parallel
    mv $tumor_bai_path ~/in/tumor_bam
    mv $normal_bai_path ~/in/normal_bam
    mv $reference_index_path ~/in/reference
    mv $reference_dict_path ~/in/reference

    docker load -i ${dockerimage_gatk_path}
  
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin/:/usr/local/bin -v /usr/local/reference/:/usr/local/reference/ -w /home/dnanexus broadinstitute/gatk:4.2.1.0 \
        /bin/bash /usr/local/bin/Mutect2.sh ${shard}.mutect.vcf.gz ${reference_path} ${tumor_bam_path} ${normal_bam_path} ${interval_file_path} ${shard}.filtered.mutect.vcf.gz

    
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
    echo "Value of tumor_sample_name: '${tumor_sample_name}'"
    echo "Value of reference: '${reference}'"
    echo "Value of reference_index: '${reference_index}'"
    
    dx download "$reference"
    dx download "$reference_index"
    dx download "$tumor_bam"
    dx download "$tumor_bai"
    

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
    
    ## concat, normalize, and extract tumor
    /usr/bin/bcftools concat --allow-overlaps --remove-duplicates -Oz --threads 8 /tmp/mutect/*.vcf.gz | /usr/bin/bcftools view -s $tumor_sample_name --threads 8 | bcftools norm -f $reference_name -m -any --threads 8 -Oz -o $output_name && /usr/bin/tabix $output_name

    vcf=$(dx upload "$output_name" --brief)
    dx-jobutil-add-output vcf --class=file "$vcf"

    vcf_index=$(dx upload "$output_name.tbi" --brief)
    dx-jobutil-add-output vcf_index --class=file "$vcf_index"

}
