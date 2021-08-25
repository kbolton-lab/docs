set -o pipefail
set -o errexit

    

export out=$1 # ${shard}.recal_data.table
export ref=$2 # ${reference_path}
export bam=$3 # ${bam_path}
export interval_file=$4 # ${interval_file_path}
export known_sites_files=$5 # "$known_sites_files"
export out_bam=$6

/gatk/gatk BaseRecalibrator --java-options "-Xmx16g" -O $out -R $ref -I $bam -L $interval_file $known_sites_files --use-original-qualities --verbosity ERROR

/gatk/gatk ApplyBQSR --java-options "-Xmx16g" -O $out_bam -I $bam --bqsr-recal-file $out -R $ref
# /gatk/gatk ApplyBQSR --java-options "-Xmx8g" -O $final_out.bqsr.bam -I $bam --bqsr-recal-file total.recal_data.table -R $ref

/usr/bin/samtools index $out_bam