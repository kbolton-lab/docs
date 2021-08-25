reference=$1
bams=$2
vcf_filename=$3
msk_out=$4
sample=$(bcftools query -l $vcf_filename)
nameroot=$5

/opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta $reference $bams --vcf $vcf_filename --output $msk_out --thread 32;
bgzip -f $msk_out && tabix -f $msk_out.gz

bcftools +fill-tags -Oz -o RD.vcf.gz $msk_out.gz -- -t "PON_RefDepth=sum(RD)"
bcftools +fill-tags -Oz -o RD_AD.vcf.gz RD.vcf.gz -- -t "PON_AltDepth=sum(AD)" && tabix RD_AD.vcf.gz

printf "##INFO=<ID=PON_RefDepth,Number=.,Type=Integer,Description=\"Total Ref_Depth for Normals\">\n##INFO=<ID=PON_AltDepth,Number=.,Type=Integer,Description=\"Total Ref and Alt Counts for Normals\">\n" > pileup.header;
printf "##INFO=<ID=PON_FISHER,Number=1,Type=Float,Description=\"P-value from Fisher's exact test with totals from PoN\">" > fisher.header;
printf "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name (with whitespace translated to underscores)\">" > sample.header;


bcftools view -H $vcf_filename | awk -v sampleID=$sample '{print $1, $2, $3, $4, $5, sampleID}' OFS='\t' > $sample.name;
bgzip -f $sample.name;
tabix $sample.name.gz -s1 -b2 -e2;
bcftools annotate --threads 32 -a $sample.name.gz -h sample.header -c CHROM,POS,-,REF,ALT,SAMPLE $vcf_filename -Oz -o $nameroot.sample.vcf.gz && tabix $nameroot.sample.vcf.gz

bcftools annotate --threads 32 -a RD_AD.vcf.gz -h pileup.header -c PON_RefDepth,PON_AltDepth $nameroot.sample.vcf.gz -Oz -o $nameroot.sample.pileup.vcf.gz;
## don't need index for VEP?
bcftools annotate --threads 32 -a RD_AD.vcf.gz -c PON_RefDepth,PON_AltDepth $nameroot.sample.vcf.gz -Oz -o $nameroot.sample.pileup.vcf.gz && tabix $nameroot.sample.pileup.vcf.gz;

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%AD]\n' $nameroot.sample.pileup.vcf.gz > $nameroot.fisher.input;


/usr/local/msk/bin/fisher.R $nameroot.fisher.input $nameroot.fisher.output
bgzip -f $nameroot.fisher.output
tabix -f -s1 -b2 -e2 $nameroot.fisher.output.gz
bcftools annotate -a $nameroot.fisher.output.gz -h fisher.header -c CHROM,POS,REF,ALT,-,-,-,-,PON_FISHER $nameroot.sample.pileup.vcf.gz -Oz -o $nameroot.fisherPON.vcf.gz && tabix $nameroot.fisherPON.vcf.gz