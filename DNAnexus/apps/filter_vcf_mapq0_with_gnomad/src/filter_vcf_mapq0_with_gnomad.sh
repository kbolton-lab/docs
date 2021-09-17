#!/bin/bash
# filter_vcf_mapq0 0.0.1

main() {

    echo "Value of vcf: '$vcf'"
    echo "Value of exclude_vcf: '$exclude_vcf'"
    echo "Value of exclude_vcf_index: '$exclude_vcf_index'"
    echo "Value of tumor_bam: '$tumor_bam'"
    echo "Value of tumor_bai: '$tumor_bai'"
    echo "Value of vcf_index: '$vcf_index'"
    echo "Value of threshold: '$threshold'"

    dx-download-all-inputs --parallel
    # mv $tumor_bai_path ~/in/tumor_bam
    # mv $vcf_index_path ~/in/vcf
    
    mv $tumor_bam_path .
    mv $vcf_path .
    mv $exclude_vcf_path .
    mv $tumor_bai_path .
    mv $vcf_index_path .
    mv $exclude_vcf_index_path .
    ## get <eid>_23153_0_0 from <eid>_23153_0_0.bqsr.bam
    eid_nameroot=$(echo $vcf_name | cut -d'.' -f1)
    echo "***********"

    # adding pre-gnomad to filter with gnomad first

    docker load -i ${dockerimage_bcf_sam_tab_path}
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/mapq0:/usr/local/mapq0 -w /home/dnanexus kboltonlab/bst \
        /bin/bash /usr/local/mapq0/mapq0_vcf_filter_samtools.sh "$vcf_name" "$tumor_bam_name" "$exclude_vcf_name" $threshold

    mapq0_vcf=$(dx upload $eid_nameroot.mapq0.soft-filtered.vcf.gz --brief)
    mapq0_vcf_index=$(dx upload $eid_nameroot.mapq0.soft-filtered.vcf.gz.tbi --brief)

    dx-jobutil-add-output mapq0_vcf "$mapq0_vcf" --class=file
    dx-jobutil-add-output mapq0_vcf_index "$mapq0_vcf_index" --class=file
}
