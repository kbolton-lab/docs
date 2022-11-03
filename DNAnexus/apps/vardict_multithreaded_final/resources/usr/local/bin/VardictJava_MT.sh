set -eou pipefail

REF="$1"
AF_THR="$2"
tumor_bam="$3"
bed="$4"
normal_bam="$5"
eid_nameroot="$6" # $eid_nameroot.vardict.vcf

tumor_sample_name=$(/usr/bin/samtools view -H $tumor_bam | /usr/bin/perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)
echo $tumor_sample_name
    
normal_sample_name=$(/usr/bin/samtools view -H $normal_bam | /usr/bin/perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)
echo $normal_sample_name

## subset bams
# tumor
tumor_bam_subset=$(basename $tumor_bam .bam).subset.bam
bedtools intersect -abam $tumor_bam -b $bed -header > $tumor_bam_subset && samtools index $tumor_bam_subset
# normal
normal_bam_subset=$(basename $normal_bam .bam).subset.bam
bedtools intersect -abam $normal_bam -b $bed -header > $normal_bam_subset && samtools index $normal_bam_subset


/opt/VarDictJava/build/install/VarDict/bin/VarDict -U -G $REF -f $AF_THR -N $tumor_sample_name -b "$tumor_bam_subset|$normal_bam_subset" -c 1 -S 2 -E 3 -g 4 $bed -th 64 | /opt/VarDictJava/build/install/VarDict/bin/testsomatic.R | /opt/VarDictJava/build/install/VarDict/bin/var2vcf_paired.pl -N "$tumor_sample_name|$normal_sample_name" -f $AF_THR > $eid_nameroot.vardict.vcf

/usr/bin/bgzip $eid_nameroot.vardict.vcf && /usr/bin/tabix $eid_nameroot.vardict.vcf.gz

# extract tumor before bcbio filter or both normal and tumor will be canidates for the filter
## bqsr changes to where Vardict NM is 0 for all????
/usr/bin/bcftools view -s ${tumor_sample_name} --threads 64 $eid_nameroot.vardict.vcf.gz | /usr/bin/bcftools norm -f $REF -m -any --threads 64 | /usr/bin/bcftools filter -e "((FMT/AF * FMT/DP < 3) && ( FMT/MQ < 55.0 || FMT/DP < 10 || FMT/QUAL < 30 ))" -m+ -s "BCBIO" --threads 64 -Oz -o $eid_nameroot.vardict.BCBIOfiltered.vcf.gz && /usr/bin/tabix $eid_nameroot.vardict.BCBIOfiltered.vcf.gz
