set -o pipefail
set -o errexit

# Usage: spark-submit run-example [options] example-class [example args]
# Options:
#   --master MASTER_URL         spark://host:port, mesos://host:port, yarn,
#                               k8s://https://host:port, or local (Default: local[*]).

export out_table=$1 # recal_data.table
export ref=$2 # ${reference_path}
export bam=$3 # ${bam_path}
export known_sites_files=$4 # "$known_sites_files"
export out_bam=$5
mkdir /tmp/spark

/gatk/gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=8" BQSRPipelineSpark \
    -R $ref -I ${bam} $known_sites_files \
    -O ${out_bam} --verbosity ERROR \
    -- --spark-runner LOCAL --spark-master local[*] \
    --conf spark.local.dir=/tmp/spark --spark-verbosity ERROR

# /gatk/gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=8" BaseRecalibratorSpark \
#     -R $ref -I ${bam} $known_sites_files \
#     -O ${out_table} --verbosity ERROR \
#     -- --spark-runner LOCAL --spark-master local[40] \
#     --conf spark.local.dir=/tmp/spark --spark-verbosity ERROR

# /gatk/gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=8" ApplyBQSRSpark \
# -R $ref -I $bam -bqsr $out_table \
# --static-quantized-quals 10 --static-quantized-quals 20 \
# --static-quantized-quals 30 -O $out_bam \
# -- --spark-runner LOCAL --spark-master local[40] \
# --conf spark.local.dir=/tmp/spark
