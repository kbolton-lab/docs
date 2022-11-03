for file in $(ls *.normalized.mutect.filtered.vcf.gz | head -n1); do
   prefix=$(basename $file .normalized.mutect.filtered.vcf.gz)
   bcftools filter -i 'FILTER="PASS"' $file | bcftools isec -C -w1 $HG38/vcf/af-only-gnomad.biallelic.above.005.leftalign.hg38.UKBB.vcf.gz
   tabix $file && bcftools isec -C -w1 $file $HG38/vcf/af-only-gnomad.biallelic.above.005.leftalign.hg38.UKBB.vcf.gz | 

   tabix ../$prefix.mutect.vep.annotated.vcf.gz
   bcftools filter -i 'FILTER="PASS"' $file -Oz -o $prefix.pass.mutect.gz && tabix -f $prefix.pass.mutect.gz && bcftools isec -C -w1 $prefix.pass.mutect.gz $HG38/vcf/af-only-gnomad.biallelic.above.005.leftalign.hg38.UKBB.vcf.gz ../$prefix.mutect.vep.annotated.vcf.gz