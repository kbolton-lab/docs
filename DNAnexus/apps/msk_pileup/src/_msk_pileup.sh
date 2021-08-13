set -eou pipefail

export sample_name="$2"
export bam_path="$3"
# optionally can get samples name from but this will be different and maybe cannot outputBinding
# export sample_name=`samtools view -H test.bam | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq`
echo "$sample_name:$bam_path"

/opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta "$1" --bam "$sample_name:$bam_path"  --vcf "$4" --output "$sample_name.pileup.vcf" --maq "$5" --baq "$6" --thread "$7";

bgzip $sample_name.pileup.vcf && tabix $sample_name.pileup.vcf.gz

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%RD]\t[%AD]\n' $sample_name.pileup.vcf.gz > $sample_name.pileup.txt

"UKB_5181257_0229406506:1006108_23153_0_0.bam"
"UKB_3716931_245493506:1010590_23153_0_0.bam"
"UKB_2797047_232997019:1043126_23153_0_0.bam"


normal_bams=( 1006108_23153_0_0.bam 1006108_23153_0_0.bam.bai 1010590_23153_0_0.bam 1010590_23153_0_0.bam.bai 1043126_23153_0_0.bam 1043126_23153_0_0.bam.bai )
vcf_name=UKB_5666276_233029566.mutect.tumor.vcf.gz
msk_out="dnanexus.pileup.vcf"

##############################
bams=()
# for bam in "${normal_bams[@]}"; do
#     dx download $bam 
#     sample_name=$(/usr/bin/samtools view -H $bam_name | /usr/bin/perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)
#      bams+=("--bam $sample_name:$bam_name")
# done
reference=/storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/GRCh38_full_analysis_set_plus_decoy_hla.fa
for bam in "${normal_bams[@]}"; do
    sample_name=$(/usr/bin/samtools view -H $bam | /usr/bin/perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)
     bams+=("--bam $sample_name:$bam")
done
echo ${bams[@]}

## if gzipped
sample=$(bcftools query -l $vcf)
filename=$(basename -- "$vcf")
extension="${filename##*.}"
if [[ $extension == "gz" ]]; then
    gunzip $vcf
    filename="${filename%.*}"
fi

## get basename
extension="${filename##*.}"
nameroot="${filename%.*}"

# /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta $reference --bam "UKB_5181257_0229406506:1006108_23153_0_0.bam" --bam "UKB_3716931_245493506:1010590_23153_0_0.bam" --bam "UKB_2797047_232997019:1043126_23153_0_0.bam" --vcf $filename --output $msk_out --maq 5 --baq 5 --thread 8;
/opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta $reference ${bams[@]} --vcf $filename --output $msk_out --thread 8;
bgzip -f $msk_out && tabix -f $msk_out.gz

bcftools +fill-tags -Oz -o RD.vcf.gz $msk_out.gz -- -t "PON_RefDepth=sum(RD)"
bcftools +fill-tags -Oz -o RD_AD.vcf.gz RD.vcf.gz -- -t "PON_AltDepth=sum(AD)" && tabix RD_AD.vcf.gz

printf "##INFO=<ID=PON_RefDepth,Number=.,Type=Integer,Description=\"Total Ref_Depth for Normals\">\n##INFO=<ID=PON_AltDepth,Number=.,Type=Integer,Description=\"Total Ref and Alt Counts for Normals\">\n" > pileup.header;
printf "##INFO=<ID=PON_FISHER,Number=1,Type=Float,Description=\"P-value from Fisher's exact test with totals from PoN\">" > fisher.header;
printf "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample name (with whitespace translated to underscores)\">" > sample.header;

sample=$(bcftools query -l $vcf)
bcftools view -H $vcf | awk -v sampleID=$sample '{print $1, $2, $3, $4, $5, sampleID}' OFS='\t' > $sample.name;
bgzip $sample.name;
tabix $sample.name.gz -s1 -b2 -e2;
bcftools annotate -a $sample.name.gz -h sample.header -c CHROM,POS,-,REF,ALT,SAMPLE $filename -Oz -o $nameroot.sample.vcf.gz && tabix $nameroot.sample.vcf.gz

bcftools annotate -a RD_AD.vcf.gz -h pileup.header -c PON_RefDepth,PON_AltDepth $nameroot.sample.vcf.gz -Oz -o $nameroot.sample.pileup.vcf.gz;
bcftools annotate -a RD_AD.vcf.gz -c PON_RefDepth,PON_AltDepth $nameroot.sample.vcf.gz -Oz -o $nameroot.sample.pileup.vcf.gz && tabix $nameroot.sample.pileup.vcf.gz;

######
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltCounts\t[%AD]\n' $name.sample.pileup.vcf > $name.fisher.input;