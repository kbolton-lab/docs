#!/bin/bash
# bqsr 0.0.1

main() {
    set -ex -o pipefail
    sudo ln -s /usr/bin/python3 /usr/bin/python

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

    ## mv everything to home
    mv $bam_path .
    mv $bam_index_path .
    mv $reference_path .
    mv $reference_index_path .
    mv $reference_dict_path .
   
    mkdir /tmp/spark

    /gatk/gatk --java-options "-Xmx32G -XX:+UseParallelGC -XX:ParallelGCThreads=16" BQSRPipelineSpark \
        -R $reference_name \
        -I $bam_name \
        $known_sites_files \
        -O $bam_prefix.bqsr.bam --verbosity ERROR \
        -- --spark-runner LOCAL --spark-master local[8] \
        --conf spark.local.dir=/tmp/spark


 
    bam_out=$(dx upload "$bam_prefix.bqsr.bam" --brief)
    dx-jobutil-add-output bam_out --class=file "$bam_out" 

    bam_out_index=$(dx upload "$bam_prefix.bqsr.bam.bai" --brief)
    dx-jobutil-add-output bam_out_index --class=file "$bam_out_index" 

}


