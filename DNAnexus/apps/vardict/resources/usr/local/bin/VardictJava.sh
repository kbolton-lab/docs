set -eou pipefail
          
REF="$1"
AF_THR="$2"
tumor_bam="$3"
tumor_sample_name="$4"
bed="$5"
normal_bam="$6"
normal_sample_name="$7"
out="$8"

/opt/VarDictJava/build/install/VarDict/bin/VarDict -G $REF -f $AF_THR -N $tumor_sample_name -b "$tumor_bam|$normal_bam" -c 1 -S 2 -E 3 -g 4 $bed | /opt/VarDictJava/build/install/VarDict/bin/testsomatic.R | /opt/VarDictJava/build/install/VarDict/bin/var2vcf_paired.pl -N "$tumor_sample_name|$normal_sample_name" -f $AF_THR > $out

bgzip $out && tabix $out.gz