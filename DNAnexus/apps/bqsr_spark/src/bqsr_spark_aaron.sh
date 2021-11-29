#!/bin/bash
# bqsr 0.0.1

main() {
    set -ex -o pipefail
    

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
    mkdir /tmp/spark
    /gatk/gatk --java-options "-Xmx124G -XX:+UseParallelGC -XX:ParallelGCThreads=32" BQSRPipelineSpark \
        -R $reference_name \
        -I $bam_name \
        $known_sites_files \
        -O $bam_prefix.bqsr.bam --verbosity ERROR \
        -- --spark-runner LOCAL --spark-master local[16] \
        --conf spark.local.dir=/tmp/spark
    # /gatk/gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=8" BQSRPipelineSpark \
    #     -R hdfs://10.0.3.103:4040/data/$reference_name \
    #     -I hdfs://10.0.3.103:4040/data/$bam_name \
    #     $known_sites_files \
    #     -O "hdfs://10.0.3.103:4040/$bam_prefix.bqsr.bam" --verbosity ERROR \
    #     -- --spark-runner SPARK --spark-master spark://23.195.26.187:7077 \
    #     --conf spark.local.dir=/tmp/spark \
    #     --executor-cores $spark_executor_cores --executor-memory $spark_executor_memory
    # /gatk/gatk --java-options "-XX:+UseParallelGC -XX:ParallelGCThreads=8" BQSRPipelineSpark \
    #     -R hdfs://10.0.3.103:4040/data/$reference_name \
    #     -I hdfs://10.0.3.103:4040/data/$bam_name \
    #     $known_sites_files \
    #     -O "hdfs://10.0.3.103:4040/$bam_prefix.bqsr.bam" --verbosity ERROR \
    #     -- --spark-runner SPARK --spark-master spark://23.195.26.187:7077 \
    #     --conf spark.local.dir=/tmp/spark \
    #     --executor-cores $spark_executor_cores --executor-memory $spark_executor_memory \
    #     --driver-memory 4g



 
    # bam_out=$(dx upload "$bam_prefix.bqsr.bam" --brief)
    # dx-jobutil-add-output bam_out --class=file "$bam_out" 

    # bam_out_index=$(dx upload "$bam_prefix.bqsr.bam.bai" --brief)
    # dx-jobutil-add-output bam_out_index --class=file "$bam_out_index" 

}


