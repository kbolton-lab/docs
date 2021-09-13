set -o pipefail
set -o errexit


export reference_name=$1 # ${reference_path}
export bam_name=$2 # ${bam_path}
export known_sites_files=$3 # "$known_sites_files"
export bam_prefix=$4

/gatk/gatk --java-options "-Xmx124G -XX:+UseParallelGC -XX:ParallelGCThreads=32" BQSRPipelineSpark \
        -R $reference_name \
        -I $bam_name \
        $known_sites_files \
        -O $bam_prefix.bqsr.bam --verbosity ERROR \
        -- --spark-runner LOCAL --spark-master local[16] \
        --conf spark.local.dir=/tmp/spark