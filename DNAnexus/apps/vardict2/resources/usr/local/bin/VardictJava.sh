set -eou pipefail
          
REF="$1"
AF_THR="$2"
tumor_bam="$3"
tumor_sample_name="$4"
bed="$5"
normal_bam="$6"
normal_sample_name="$7"
out="$8"
interval_file_name="$9"

## subset bams
# tumor
# tumor_bam_subset=$(basename $tumor_bam .bam).subset.bam
# bedtools intersect -abam $tumor_bam -b $bed -header > $tumor_bam_subset && samtools index $tumor_bam_subset
# normal
# normal_bam_subset=$(basename $normal_bam .bam).subset.bam
# bedtools intersect -abam $normal_bam -b $bed -header > $normal_bam_subset && samtools index $normal_bam_subset
bedtools makewindows -b $bed -w 1150 -s 1000 > hg38.interval.1K.bed

# JAVA_OPTS=-Xmx400g
# /opt/VarDictJava/build/install/VarDict/bin/VarDict -fisher -G $REF -f $AF_THR -N $tumor_sample_name -b "$tumor_bam|$normal_bam" -c 1 -S 2 -E 3 -g 4 hg38.interval.50K.bed -th 32 | /opt/VarDictJava/build/install/VarDict/bin/testsomatic.R | /opt/VarDictJava/build/install/VarDict/bin/var2vcf_paired.pl -N "$tumor_sample_name|$normal_sample_name" -f $AF_THR > $out
/opt/VarDictJava/build/install/VarDict/bin/VarDict -fisher -G $REF -f $AF_THR -N $tumor_sample_name -b "$tumor_bam|$normal_bam" -c 1 -S 2 -E 3 -g 4 hg38.interval.1K.bed -th 16 > $interval_file_name.txt
cat $interval_file_name.txt | /opt/VarDictJava/build/install/VarDict/bin/var2vcf_paired.pl -N "$tumor_sample_name|$normal_sample_name" -f $AF_THR > $out

bgzip $out && tabix $out.gz