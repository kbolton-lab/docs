#!/bin/bash

set -ex -o pipefail
echo "brian"

vcf=$1 # pre-gnomad vcf
bam=$2
exclude_vcf=$3
mapq0perc=$4
eid_nameroot=$5

# adding pre-gnomad to filter with gnomad first
# grab sites that don't already have the MQ0 field
bcftools isec -C -w1 ${vcf} ${exclude_vcf} -Oz -o gnomad.filtered.vcf.gz && tabix gnomad.filtered.vcf.gz

zgrep -v "^#" gnomad.filtered.vcf.gz | grep -v "MQ0" | cut -f 1,2 | while read chr pos; do
    /usr/local/bin/samtools stats -d -@8 "$bam" $chr:$pos-$pos > stats
    mapq0=$(grep "reads MQ0:" stats | cut -f3); printf "$chr\t$pos\t$mapq0\t" >> mapq0counts
    # sequences/reads
    grep -P "SN\tsequences:" stats | cut -f3 >> mapq0counts
done

printf "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Number of MAPQ == 0 reads covering this record\">\n##INFO=<ID=samtools_DP,Number=1,Type=Integer,Description=\"Samtools depth at this position\">\n" > MQ0.header;

bgzip -f mapq0counts
tabix mapq0counts.gz -s1 -b2 -e2;
bcftools annotate --threads 8 -a mapq0counts.gz -h MQ0.header -c CHROM,POS,MQ0,samtools_DP gnomad.filtered.vcf.gz | bcftools filter -m+ -e "((INFO/MQ0) / (INFO/samtools_DP)) > $mapq0perc" -s "MQ0" --threads 8 -Oz -o $eid_nameroot.mapq0.soft-filtered.vcf.gz && tabix $eid_nameroot.mapq0.soft-filtered.vcf.gz
