bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltCounts\t[%AD]\n' $vcf_path > $vcf_prefix.fisher.input;


LC_ALL=C.UTF-8 Rscript --vanilla /usr/local/R/bin/fisher.R $vcf_prefix.fisher.input $vcf_prefix.fisher.output
bgzip -f $vcf_prefix.fisher.output
tabix -f -s1 -b2 -e2 $vcf_prefix.fisher.output.gz
bcftools annotate -a $vcf_prefix.fisher.output.gz -h fisher.header -c CHROM,POS,REF,ALT,-,-,-,-,PON_FISHER $vcf_prefix.sample.pileup.vcf -Oz -o $vcf_prefix.fisherPON.vcf.gz && tabix $vcf_prefix.fisherPON.vcf.gz


## just include p-value instead of hard filter
#bcftools filter -i "INFO/PON_FISHER<$p_value" $name.pileup.fisherPON.vcf.gz -Oz -o $name.filtered.pileup.fisherPON.vcf.gz