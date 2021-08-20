set -o pipefail
set -o errexit

out=$1
inputs=$2
bam=$3
final_out=$4
ref=$5

/gatk/gatk GatherBQSRReports --java-options "-Xmx8g" -O $out $inputs

/gatk/gatk ApplyBQSR --java-options "-Xmx8g" -O $final_out -I $bam --bqsr-recal-file $out -R $ref
