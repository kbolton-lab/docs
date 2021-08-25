#!/bin/bash
# bqsr 0.0.1

main() {
    set -ex -o pipefail
    # locate start-master.sh
    which python3
    # ls /sbin
    # locate start-master.sh
    sudo ln -s /usr/bin/python3 /usr/bin/python
    # echo $HADOOP_HOME # /cluster/hadoop
    # echo $SPARK_HOME # /cluster/spark
    # which hadoop # /cluster/hadoop/bin/hadoop
    # which spark-submit # /cluster/spark/bin/spark-submit
    # ls /cluster
    # $HADOOP_HOME/bin/hadoop -h
    # $SPARK_HOME/bin/spark-submit -h
    
    #exit
    
    # echo "Value of bam: '${bam}'"
    # echo "Value of bam_index: '${bam_index}'"
    # echo "Value of known_sites: '${known_sites[@]}'"
    # echo "Value of known_sites_index: '${known_sites_index[@]}'"
    $HADOOP_HOME/bin/hadoop fs -mkdir -p /data
    dx-download-all-inputs --parallel

    known_sites_files=""
    for file in "${known_sites_path[@]}"; do
        known_sites_filename=$(basename -- "$file")
        ## mv all known sites files to home
        mv $file .
        known_sites_files="$known_sites_files --known-sites $known_sites_filename" 
        # $HADOOP_HOME/bin/hadoop fs -put $file /data/$known_sites_filename
        # known_sites_files="$known_sites_files --known-sites hdfs://10.0.3.103:4040/data/$known_sites_filename"   
    done

    ## the index all have different paths, don't need file names just move to where their file is locates which is /home/dnanexus
    for file in "${known_sites_index_path[@]}"; do
        ## mv all indexes to home
        mv $file .
        # $HADOOP_HOME/bin/hadoop fs -put $file /data
    done

    ## mv everything to hdfs
    mv $bam_path .
    mv $bam_index_path .
    mv $reference_path .
    mv $reference_index_path .
    mv $reference_dict_path .
   
    # $HADOOP_HOME/bin/hadoop fs -put $bam_path /data/$bam_name
    # $HADOOP_HOME/bin/hadoop fs -put $bam_index_path /data/$bam_index_name
    # $HADOOP_HOME/bin/hadoop fs -put $reference_path /data/$reference_name
    # $HADOOP_HOME/bin/hadoop fs -put $reference_index_path /data/$reference_index_name
    # $HADOOP_HOME/bin/hadoop fs -put $reference_dict_path /data/$reference_dict_name
    
    #mv $dockerimage_gatk_path .
    #dx download "${dockerimage_gatk}"
    #docker load -i ${dockerimage_gatk_path}
    hdfs dfs -ls /data
    hadoop fs -ls /data

    spark_executor_memory=124g # Could be parameterized as an app input
    spark_executor_cores=16  # Could be parameterized as an app input
    mkdir /tmp/spark

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
    /gatk/gatk --java-options "-Xmx124G -XX:+UseParallelGC -XX:ParallelGCThreads=32" BQSRPipelineSpark \
        -R $reference_name \
        -I $bam_name \
        $known_sites_files \
        -O $bam_prefix.bqsr.bam --verbosity ERROR \
        -- --spark-runner LOCAL --spark-master local[16] \
        --conf spark.local.dir=/tmp/spark


    # docker run --rm \
    # -v /home/dnanexus:/home/dnanexus \
    # -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 \
    # -v /usr/local/helper_script/bin:/usr/local/helper_script/bin \
    # -v $HADOOP_HOME:$HADOOP_HOME \
    # --env HADOOP_HOME=$HADOOP_HOME \
    # -v $SPARK_HOME:$SPARK_HOME \
    # --env SPARK_HOME=$SPARK_HOME \
    # -v /cluster:/cluster \
    # -w /home/dnanexus \
    # broadinstitute/gatk:4.2.1.0 /bin/bash /usr/local/helper_script/bin/bqsr_spark_helper2.sh \
    #     $spark_executor_memory $spark_executor_cores /hdfs/data/$reference_index_name

    ## this errored out on job-G4Vyyy0J6XG97b1xFJf2BvPZ
    # /cluster/log_collector.sh /home/dnanexus/out/cluster_runtime_logs_tarball
    cat /cluster/log_collector.sh
    ls 
    ls /cluster
 
    bam_out=$(dx upload "$bam_prefix.bqsr.bam" --brief)
    dx-jobutil-add-output bam_out --class=file "$bam_out" 

    bam_out_index=$(dx upload "$bam_prefix.bqsr.bam.bai" --brief)
    dx-jobutil-add-output bam_out_index --class=file "$bam_out_index" 

}


