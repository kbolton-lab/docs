set -o pipefail
set -o errexit

export out=$1
export ref=$2
export bam=$3
export interval_file=$4
export known_sites_files=$5

/gatk/gatk BaseRecalibrator --java-options "-Xmx8g" -O $out -R $ref -I $bam -L $interval_file $known_sites_files --use-original-qualities --verbosity ERROR
