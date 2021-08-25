set -o pipefail
set -o errexit

mkdir /tmp/spark
spark_executor_memory=$1
spark_executor_cores=$2
$HADOOP_HOME/bin/hadoop -h
$SPARK_HOME/bin/spark-submit -h

ref_index=$3
export PATH=$PATH:$HADOOP_HOME/bin/:$SPARK_HOME/bin/
$HADOOP_HOME/bin/hadoop fs -mkdir -p /hdfs/data
$HADOOP_HOME/bin/hadoop fs -put $ref_index /hdfs/data/$ref_index
which hadoop
which spark-submit
ls /hdfs/data/
ls /hdfs/data/$ref_index


/gatk/gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=8" BQSRPipelineSpark \
    -R $ref -I ${bam} $known_sites_files \
    -O ${out_bam} --verbosity ERROR \
    -- --spark-runner LOCAL --spark-master local[*] \
    --conf spark.local.dir=/tmp/spark --spark-verbosity ERROR


