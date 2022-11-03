#!/bin/bash
# mutect 0.0.1

main() {
    set -ex -o pipefail
    echo ${interval_list_path}
    echo "BRIAN"
    dx-download-all-inputs --parallel
    mv $tumor_bai_path ~/in/tumor_bam
    mv $reference_index_path ~/in/reference
    mv $reference_dict_path ~/in/reference

    docker load -i ${dockerimage_gatk_path}
    ## get <eid>_23153_0_0 from <eid>_23153_0_0.bqsr.bam
    eid_nameroot=$(echo $tumor_bam_name | cut -d'.' -f1)
    output_name="$eid_nameroot.haplotype_caller.vcf.gz"
  
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin/:/usr/local/bin -v /usr/local/reference/:/usr/local/reference/ -w /home/dnanexus broadinstitute/gatk:4.2.1.0 \
        /bin/bash /usr/local/bin/ht_call.sh ${reference_path} ${tumor_bam_path} $output_name ${interval_list_path} 

 
    output_name_norm="$eid_nameroot.normalized.haplotype_caller.vcf.gz"
    
    ## normalize
    /usr/bin/bcftools norm -f $reference_path -m -any --threads 8 -Oz -o $output_name_norm $output_name && /usr/bin/tabix $output_name_norm

    /usr/local/bin/bcftools norm -f $reference_path -m -any --threads 8 1199870_23153_0_0.raw.vcf.gz | ZL

    vcf=$(dx upload "$output_name_norm" --brief)
    dx-jobutil-add-output vcf --class=file "$vcf"

    vcf_index=$(dx upload "$output_name_norm.tbi" --brief)
    dx-jobutil-add-output vcf_index --class=file "$vcf_index"

}
