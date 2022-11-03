#!/bin/bash

set -ex -o pipefail
echo "brian"

vcf=$1
bam=$2
mapq0perc=$3
eid_nameroot=$(echo $vcf | cut -d'.' -f1)

#grab sites that don't already have the MQ0 field
zgrep -v "^#" "$vcf" | grep -v "MQ0" | cut -f 1,2 | while read chr pos; do
    mapq0=$(samtools stats -d -@8 "$bam" $chr:$pos-$pos | grep "reads MQ0:" | cut -f3); printf "$chr\t$pos\t$mapq0\n" >> mapq0counts
done

printf "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Number of MAPQ == 0 reads covering this record\">" > MQ0.header;

bgzip -f mapq0counts
tabix mapq0counts.gz -s1 -b2 -e2;
bcftools annotate --threads 8 -a mapq0counts.gz -h MQ0.header -c CHROM,POS,MQ0 $vcf | bcftools filter -m+ -e "((INFO/MQ0) / (FMT/DP)) > $mapq0perc" -s "MQ0" --threads 8 -Oz -o $eid_nameroot.mapq0.soft-filtered.gz && tabix $eid_nameroot.mapq0.soft-filtered.gz