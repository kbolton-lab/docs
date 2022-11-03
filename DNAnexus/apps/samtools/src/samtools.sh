#!/bin/bash
# cram2bam 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://documentation.dnanexus.com/developer for tutorials on how
# to modify this file.

main() {
    set -ex -o pipefail
    echo "Value of cram_in: '$cram_in'"
    echo "Value of crai_in_path: '$crai_in_path'"
    echo "Value of ref_fasta: '$ref_fasta'"

    dx-download-all-inputs --parallel
    mv $crai_in_path /home/dnanexus/in/cram_in
    bam_out="${cram_in_prefix}.bam"
	

    #/usr/bin/samtools_helper.sh $ref_fasta_path $cram_in_path $bam_out $region
    /usr/bin/samtools view -@16 -b -T $ref_fasta_path $cram_in_path -o $bam_out
    /usr/bin/samtools index $bam_out;
    
    # save fastq reads in separate R1 and R2 files

    /usr/bin/samtools fastq -@16 $bam_out \
        -1 ${cram_in_prefix}_R1.fastq.gz \
        -2 ${cram_in_prefix}_R2.fastq.gz \
        -0 /dev/null -s /dev/null -n

    # dx run applet-GFp8G80JP6z1jVP78yGQB4qZ -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/samtools/inputs.json --destination $KELLY:/Brian/test/ -y
    # dx run pbgermline -ireads_fastqgzs=file-GFp6fV8JF2JkQZqj78J6gVyX -ireads2_fastqgzs=file-GFp6fX0JF2JjVGjQBvGPkv4j -igenome_fastagz=file-G3jKk3QJ6XG9KY8680qg9j8b -igenomeindex_targz=file-GFp96p0J6XG998yGJ5f88x58 -ilicense_file_bin=file-GFpB1jjJQ28KZX9G4vvQQYPf --destination $KELLY:/Brian/test/ -y --extra-args '{"regionalOptions": { "aws:us-west-2": { "systemRequirements": { "*": { "instanceType": "mem2_ssd1_gpu1_x32" }}}}}'
    # bai_out=$(dx upload ${bam_out}.bai --brief)
    # bam_out=$(dx upload $bam_out --brief)
   
    # dx-jobutil-add-output bam_out "$bam_out" --class=file
    # dx-jobutil-add-output bai_out "$bai_out" --class=file
    fq1=$(dx upload ${cram_in_prefix}_R1.fastq.gz --brief)
    fq2=$(dx upload ${cram_in_prefix}_R2.fastq.gz --brief)
   
    dx-jobutil-add-output fq1 "$fq1" --class=file
    dx-jobutil-add-output fq2 "$fq2" --class=file
}