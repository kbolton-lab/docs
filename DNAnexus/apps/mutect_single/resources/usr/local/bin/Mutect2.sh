set -o pipefail
set -o errexit

export tumor_bam="$3"
export normal_bam="$4"

NORMAL=$(samtools view -H $normal_bam | /usr/bin/perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)
TUMOR=$(samtools view -H $tumor_bam | /usr/bin/perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)

/gatk/gatk Mutect2 --java-options "-Xmx8g" -O $1 -R $2 -I $3 -tumor "$TUMOR" -I $4 -normal "$NORMAL" -L $5 #Running Mutect2.

/gatk/gatk FilterMutectCalls --java-options "-Xmx8g" -R $2 -V $1 -O $6 #Running FilterMutectCalls on the output vcf.
