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


java -jar /varscan/VarScan.v2.4.2.jar somatic \
    <(/usr/bin/samtools mpileup --no-baq -l ${bed} -f ${REF} ${normal_bam} ${tumor_bam}) "${out}" \
    --strand-filter 0 \
    --min-coverage 8 \
    --min-var-freq 0.005 \
    --p-value 0.99 \
    --mpileup 1 \
    --output-vcf


/usr/bin/bgzip ${out}.snp.vcf && /usr/bin/tabix ${out}.snp.vcf.gz
/usr/bin/bgzip ${out}.indel.vcf && /usr/bin/tabix ${out}.indel.vcf.gz
/usr/bin/bcftools concat -a -D ${out}.snp.vcf.gz ${out}.indel.vcf.gz -Oz -o ${out}.concat.vcf.gz && /usr/bin/tabix ${out}.concat.vcf.gz

printf "TUMOR ${tumor_sample_name}\nNORMAL ${normal_sample_name}\n" > sample_update.txt
/usr/bin/bcftools reheader ${out}.concat.vcf.gz -s sample_update.txt -o ${out}.concat.reheader.vcf.gz && /usr/bin/tabix -p vcf ${out}.concat.reheader.vcf.gz