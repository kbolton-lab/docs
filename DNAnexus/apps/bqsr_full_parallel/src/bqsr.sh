#!/bin/bash
# bqsr 0.0.1

main() {
    set -ex -o pipefail
    echo "Value of bam: '${bam}'"
    echo "Value of bam_index: '${bam_index}'"
    echo "Value of known_sites: '${known_sites[@]}'"
    echo "Value of known_sites_index: '${known_sites_index[@]}'"
    echo "Value of interval_lists_array: '${interval_lists_array[@]}'"

   
    interval_files=()
    for file in "${interval_lists_array[@]}"; do
        interval_files+=("-iarray_of_scattered_input=${file}")
    done

    known_sites_files=()
    for file in "${known_sites[@]}"; do
        known_sites_files+=("-iknown_sites=${file}")
    done

    known_sites_index_files=()
    for file in "${known_sites_index[@]}"; do
        known_sites_index_files+=("-iknown_sites_index=${file}")
    done

    map_job=$(dx-jobutil-new-job map \
        "${interval_files[@]}" \
        -i "reference=${reference}" \
        -i "reference_index=${reference_index}" \
        -i "reference_dict=${reference_dict}" \
        -i "bam=${bam}" \
        -i "bam_index=${bam_index}" \
        "${known_sites_files[@]}" \
        "${known_sites_index_files[@]}" \
        -i "dockerimage_gatk=${dockerimage_gatk}")

    
    #output_name="$bam_prefix.bqsr.bam"
    postprocess_job=$(dx-jobutil-new-job postprocess \
        -i process_outputs:array:jobref=${map_job}:process_outputs \
        -i "bam=${bam}" \
        -i "bam_index=${bam_index}" \
        -i "reference=${reference}" \
        -i "reference_index=${reference_index}" \
        -i "reference_dict=${reference_dict}" \
        -i "dockerimage_gatk=${dockerimage_gatk}" \
        --depends-on $map_job)

    dx-jobutil-add-output bam_out --class=jobref "$postprocess_job":bam_out 
    dx-jobutil-add-output bam_out_index --class=jobref "$postprocess_job":bam_out_index
}

map() {
    set -ex -o pipefail
    
    echo "Value of array_of_scattered_input: '${array_of_scattered_input[@]}'" 
    echo "Value of known_sites: '${known_sites[@]}'"
    echo "Value of known_sites_index: '${known_sites_index[@]}'"

    known_sites_files=()
    for file in "${known_sites[@]}"; do
        known_sites_files+=("-iknown_sites=${file}")
    done

    known_sites_index_files=()
    for file in "${known_sites_index[@]}"; do
        known_sites_index_files+=("-iknown_sites_index=${file}")
    done

    process_jobs=()
    i=1
    for interval_file in "${array_of_scattered_input[@]}"
    do
        echo "Value of interval_file: '${interval_file}'"
        
        process_jobs+=($(dx-jobutil-new-job process -i "interval_file=${interval_file}" -i "reference=${reference}" -i "reference_index=${reference_index}" -i "reference_dict=${reference_dict}" -i "bam=${bam}" -i "bam_index=${bam_index}" -i shard:string=shard-$i  "${known_sites_files[@]}"  "${known_sites_index_files[@]}" -i "dockerimage_gatk=${dockerimage_gatk}"))
        echo "${process_jobs[@]}"
        i=$((i+1))
        echo "$i.after"
    done

    for process_job in "${process_jobs[@]}"
    do
        dx-jobutil-add-output process_outputs --class=array:jobref "${process_job}:table" --array
    done

}

process() {
    set -ex -o pipefail

    echo "Value of known_sites: '${known_sites[@]}'"
    echo "Value of known_sites_index: '${known_sites_index[@]}'"
    
    dx-download-all-inputs --parallel

    known_sites_files=""
    for file in "${known_sites_path[@]}"; do
        known_sites_filename=$(basename -- "$file")
        mv $file .
        known_sites_files="$known_sites_files --known-sites $known_sites_filename"   
    done

    ## the index all have different paths, don't need file names just move to where their file is locates which is /home/dnanexus
    for file in "${known_sites_index_path[@]}"; do
        mv $file .
    done

    mv $bam_index_path ~/in/bam
    mv $reference_index_path ~/in/reference
    mv $reference_dict_path ~/in/reference
    docker load -i ${dockerimage_gatk_path}
  
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin/:/usr/local/bin -w /home/dnanexus broadinstitute/gatk:4.2.1.0 \
        /bin/bash /usr/local/bin/run_bqsr.sh ${shard}.recal_data.table ${reference_path} ${bam_path} ${interval_file_path} "$known_sites_files"
 

    table=$(dx upload ${shard}.recal_data.table --brief)
    dx-jobutil-add-output table --class=file "$table" 

}

postprocess() {
    set -ex -o pipefail

    echo "Value of bam: '${bam}'"
    echo "Value of bam_index: '${bam_index}'"
    echo "Value of process_outputs: '${process_outputs[@]}'"
    echo "Value of output_name: '${output_name}'"

    dx download "$bam"
    dx download "$bam_index"
    dx download "$reference"
    dx download "$reference_index"
    dx download "$reference_dict" 
    dx download "$dockerimage_gatk"

    # download vcfs and vcf indices
    mkdir -p /tmp/bqsr_tables
    cd /tmp/bqsr_tables
    for process_output in "${process_outputs[@]}" 
    do
        dx download "$process_output"
        echo "$process_output_path"
    done
    ls -1sh /tmp/bqsr_tables/*.recal_data.table

    inputs=""
    for file in $(ls /tmp/bqsr_tables/*.recal_data.table); do
        inputs="$inputs -I $file"   
    done

    cd /home/dnanexus
    output_name="$bam_prefix.bqsr.bam"
    echo $output_name
    final_out=$(basename $bam_name .bam)



    docker load -i ${dockerimage_gatk_name}
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin/:/usr/local/bin -v /tmp/bqsr_tables:/tmp/bqsr_tables -w /home/dnanexus broadinstitute/gatk:4.2.1.0 /bin/bash /usr/local/bin/gather_bqsr.sh "$inputs" $bam_name $final_out $reference_name

    bam_out=$(dx upload "$final_out.bqsr.bam" --brief)
    dx-jobutil-add-output bam_out --class=file "$bam_out" 

    bam_out_index=$(dx upload "$final_out.bqsr.bam.bai" --brief)
    dx-jobutil-add-output bam_out_index --class=file "$bam_out_index"

}
