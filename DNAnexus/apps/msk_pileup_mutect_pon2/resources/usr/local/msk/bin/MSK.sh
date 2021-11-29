set -ex -o pipefail

reference=$1
bams=$2
vcf_filename=$3
msk_out=$4
sample=$(bcftools query -l $vcf_filename)
nameroot=$5
vcf2PON=$6
caller=$7
p_value=$8

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
# bcftools annotate --threads 32 -a RD_AD.vcf.gz -c PON_RefDepth,PON_AltDepth $nameroot.sample.vcf.gz -Oz -o $nameroot.sample.pileup.vcf.gz && tabix $nameroot.sample.pileup.vcf.gz;

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%AD]\n' $nameroot.sample.pileup.vcf.gz > $nameroot.fisher.input;


/usr/local/msk/bin/fisher.R $nameroot.fisher.input $nameroot.fisher.output
bgzip -f $nameroot.fisher.output
tabix -f -s1 -b2 -e2 $nameroot.fisher.output.gz
bcftools annotate -a $nameroot.fisher.output.gz -h fisher.header -c CHROM,POS,REF,ALT,-,-,-,-,PON_FISHER $nameroot.sample.pileup.vcf.gz -Oz -o $nameroot.fisherPON.vcf.gz && tabix $nameroot.fisherPON.vcf.gz

export name="$nameroot"."$caller".msk.pon2.vcf.gz
printf "##INFO=<ID=PON_2AT2_percent,Number=1,Type=Integer,Description=\"If 2 PoN samples have variant at >=2 percent\">\n" > pon2.header;
bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t1\n" $vcf2PON > normal2.txt
bgzip -f normal2.txt
tabix -f -s1 -b2 -e2 normal2.txt.gz
## now hard filtering PoN 
bcftools annotate --threads 4 -a normal2.txt.gz -h pon2.header -c CHROM,POS,REF,ALT,PON_2AT2_percent $nameroot.fisherPON.vcf.gz | bcftools filter -i "PON_FISHER <= $p_value" -Oz -o $name
tabix $name