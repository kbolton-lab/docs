#!/bin/bash
# filter_vcf_mapq0 0.0.1

main() {

    echo "Value of vcf: '$vcf'"
    echo "Value of tumor_bam: '$tumor_bam'"
    echo "Value of tumor_bai: '$tumor_bai'"
    echo "Value of vcf_index: '$vcf_index'"
    echo "Value of threshold: '$threshold'"

    dx-download-all-inputs --parallel
    # mv $tumor_bai_path ~/in/tumor_bam
    # mv $vcf_index_path ~/in/vcf
    
    mv $tumor_bam_path .
    mv $vcf_path .
    mv $tumor_bai_path .
    mv $vcf_index_path .
    ## get <eid>_23153_0_0 from <eid>_23153_0_0.bqsr.bam
    eid_nameroot=$(echo $vcf_name | cut -d'.' -f1)
    echo "***********"

    # #grab sites that don't already have the MQ0 field
    # zgrep -v "^#" "$vcf_path" | grep -v "MQ0" | cut -f 1,2 | while read chr pos; do
    #     mapq0=$(samtools stats -d -@8 "$tumor_bam_path" $chr:$pos-$pos |  grep "reads MQ0:" | cut -f3); printf "$chr\t$pos\t$mapq0\n" >> mapq0counts
    # done

    # printf "##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Number of MAPQ == 0 reads covering this record\">" > MQ0.header;

    # bgzip -f mapq0counts
    # tabix mapq0counts.gz -s1 -b2 -e2;
    # bcftools annotate --threads 8 -a mapq0counts.gz -h MQ0.header -c CHROM,POS,MQ0 $vcf_path | bcftools filter -e "((INFO/MQ0) / (FMT/DP)) > $threshold" -s "MQ0" --threads 8 -Oz -o $eid_nameroot.mapq0.soft-filtered.gz && tabix $eid_nameroot.mapq0.soft-filtered.gz

    docker load -i ${dockerimage_bcf_sam_tab_path}
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/mapq0:/usr/local/mapq0 -w /home/dnanexus kboltonlab/bst \
        /bin/bash /usr/local/mapq0/mapq0_vcf_filter_samtools.sh "$vcf_name" "$tumor_bam_name" $threshold

    mapq0_vcf=$(dx upload $eid_nameroot.mapq0.soft-filtered.gz --brief)
    mapq0_vcf_index=$(dx upload $eid_nameroot.mapq0.soft-filtered.gz.tbi --brief)

    dx-jobutil-add-output mapq0_vcf "$mapq0_vcf" --class=file
    dx-jobutil-add-output mapq0_vcf_index "$mapq0_vcf_index" --class=file
}
