#!/bin/bash

ref_fasta_path=$1
cram_in_path=$2
bam_out=$3

/usr/local/bin/samtools view -b -T $ref_fasta_path $cram_in_path -o $bam_out;
/usr/local/bin/samtools index $bam_out;
sample_name=$(samtools view -H $bam_out | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq)
