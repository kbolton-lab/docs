#!/bin/bash
# annotate_CH_pd 0.0.1


main() {
    set -ex -o pipefail
    echo "8/20/2021"

    echo "Value of mutect_vcf: '$mutect_vcf'"
    echo "Value of vardict_vcf: '$vardict_vcf'"
    echo "Value of tumor_sample_name: '$tumor_sample_name'"
    echo "Value of impact_annotation: '$impact_annotation'"
    echo "Value of intervals: '$intervals'"
    echo "Value of dockerimage_annotate_PD: '$dockerimage_annotate_PD'"

    
    dx-download-all-inputs --parallel

    if [ -z "$tumor_sample_name" ]
    then
        echo "\$tumor_sample_name is empty"
        tumor_sample_name=$(/usr/bin/bcftools query -l $mutect_vcf_path)
        echo $tumor_sample_name
    fi

    ## get <eid>_23153_0_0 from <eid>_23153_0_0.<suffix>
    eid_nameroot=$(echo $mutect_vcf_name | cut -d'.' -f1)

    ls $mutect_vcf_path
    ls $vardict_vcf_path
    ls $impact_annotation_path
    ls $topmed_annotation_path
    ls $cosmic_annotation_path
    ls $tsg_annotation_path
    ls $oncoKB_annotation_path
    ls $pd_table_annotation_path
    ls $panmyeloid_annotation_path
    ls $blacklist_annotation_path
    ls $segemental_duplications_annotation_path
    ls $simple_repeats_annotation_path
    ls $repeat_masker_annotation_path

    # if [ -n "$impact_annotation" ]
    # then
    #     dx download "$impact_annotation" -o impact_annotation
    # fi
    target_length=$(grep -v "^@" $intervals_path | awk '{print $1,$2,$3}' | awk '$4=$3-$2 {sum+=$4}END {print sum}')

    docker load -i $dockerimage_annotate_PD_path
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/annotate/:/usr/local/annotate -w /home/dnanexus kboltonlab/annotate_wes_ch:3.0 /usr/local/annotate/annotate_CH_pd_docker3_Mutect_Vardict_ponChange.R \
        -v $mutect_vcf_path,$vardict_vcf_path \
        -i $impact_annotation_path \
        -b $topmed_annotation_path \
        -c $cosmic_annotation_path \
        -T $tsg_annotation_path \
        --oncoKB-curated $oncoKB_annotation_path \
        -p $pd_table_annotation_path \
        --pan-myeloid $panmyeloid_annotation_path \
        --blacklist $blacklist_annotation_path \
        --segemental-duplications $segemental_duplications_annotation_path \
        --simple-repeats $simple_repeats_annotation_path \
        --repeat-masker $repeat_masker_annotation_path \
        --target-length $target_length \
        -e $eid_nameroot
    
    final_tsv=$(dx upload $eid_nameroot.final.tsv --brief)
    column_check=$(dx upload $eid_nameroot.columns.txt --brief)

    dx-jobutil-add-output final_tsv "$final_tsv" --class=file
    dx-jobutil-add-output column_check "$column_check" --class=file
}
