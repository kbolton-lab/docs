set -o pipefail
set -o errexit

export tumor_bam="$1"
which perl
which samtools

TUMOR=$(samtools view -H $tumor_bam | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)
echo $TUMOR