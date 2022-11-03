#!/bin/bash

ref_fasta_path=$1
cram_in_path=$2
bam_out=$3

# /usr/local/bin/samtools view -@4 -b -T $ref_fasta_path $cram_in_path -o $bam_out;
# /usr/local/bin/samtools index $bam_out;

/usr/bin/samtools view -@4 -b -T $ref_fasta_path $cram_in_path -o $bam_out;
/usr/bin/samtools index $bam_out;
sample_name=$(/usr/bin/samtools view -H $bam_out | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)
