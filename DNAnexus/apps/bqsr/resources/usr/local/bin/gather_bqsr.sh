set -o pipefail
set -o errexit

out=$1
inputs=$1
bam=$2
final_out=$3
ref=$4

echo $inputs

/gatk/gatk GatherBQSRReports --java-options "-Xmx8g" -O total.recal_data.table $inputs

/gatk/gatk ApplyBQSR --java-options "-Xmx8g" -O $final_out.bqsr.bam -I $bam --bqsr-recal-file total.recal_data.table -R $ref

samtools index $final_out.bqsr.bam
