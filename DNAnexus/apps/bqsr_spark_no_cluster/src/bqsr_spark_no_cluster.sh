#!/bin/bash
# bqsr 0.0.1

main() {
    set -ex -o pipefail

    
    dx-download-all-inputs --parallel

    known_sites_files=""
    for file in "${known_sites_path[@]}"; do
        known_sites_filename=$(basename -- "$file")
        ## mv all known sites files to home
        mv $file .
        known_sites_files="$known_sites_files --known-sites $known_sites_filename" 
    done

    ## the index all have different paths, don't need file names just move to where their file is locates which is /home/dnanexus
    for file in "${known_sites_index_path[@]}"; do
        ## mv all indexes to home
        mv $file .
    done

    ## mv everything to hdfs
    mv $bam_path .
    mv $bam_index_path .
    mv $reference_path .
    mv $reference_index_path .
    mv $reference_dict_path .
 
    docker load -i ${dockerimage_gatk_path}
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/helper_script/bin:/usr/local/helper_script/bin -w /home/dnanexus broadinstitute/gatk:4.2.1.0 \
        /bin/bash /usr/local/helper_script/bin/bqsr_spark_helper.sh $reference_name $bam_name "$known_sites_files" $bam_prefix
 
    bam_out=$(dx upload "$bam_prefix.bqsr.bam" --brief)
    dx-jobutil-add-output bam_out --class=file "$bam_out" 

    bam_out_index=$(dx upload "$bam_prefix.bqsr.bam.bai" --brief)
    dx-jobutil-add-output bam_out_index --class=file "$bam_out_index" 

}


