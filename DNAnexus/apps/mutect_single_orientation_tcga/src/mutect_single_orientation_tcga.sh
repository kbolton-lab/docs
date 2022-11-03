#!/bin/bash
# mutect 0.0.1

main() {
    set -ex -o pipefail
    echo ${interval_list_path}
    echo "BRIAN"
    dx-download-all-inputs --parallel
    mv $tumor_bai_path ~/in/tumor_bam
    mv $normal_bai_path ~/in/normal_bam
    mv $reference_index_path ~/in/reference
    mv $reference_dict_path ~/in/reference

    docker load -i ${dockerimage_gatk_path}
    ## TCGA-blah-blah-.bqsr.bam -> TCGA-blah-blah-.mutect.filtered.vcf.gz
    eid_nameroot=$tumor_bam_prefix
    output_name="$eid_nameroot.mutect.filtered.vcf.gz"
  
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin/:/usr/local/bin -v /usr/local/reference/:/usr/local/reference/ -w /home/dnanexus broadinstitute/gatk:4.2.1.0 \
        /bin/bash /usr/local/bin/Mutect2.sh "$eid_nameroot.mutect.vcf.gz" ${reference_path} ${tumor_bam_path} ${normal_bam_path} ${interval_list_path} $output_name

    tumor_sample_name=$(/usr/bin/samtools view -H ${tumor_bam_path} | /usr/bin/perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)
 
    output_name_norm="$eid_nameroot.normalized.mutect.filtered.vcf.gz"
    
    ## extract "tumor"/sample and normalize
    /usr/bin/bcftools view -s $tumor_sample_name --threads 8 $output_name | bcftools norm -f $reference_path -m -any --threads 8 -Oz -o $output_name_norm && /usr/bin/tabix $output_name_norm

    vcf=$(dx upload "$output_name_norm" --brief)
    dx-jobutil-add-output vcf --class=file "$vcf"

    vcf_index=$(dx upload "$output_name_norm.tbi" --brief)
    dx-jobutil-add-output vcf_index --class=file "$vcf_index"

}
