#!/bin/bash
# annotate_CH_pd 0.0.1


main() {
    set -ex -o pipefail
    echo "10/15/2021"
    ls /usr/local/annotate/annotate_CH_pd_docker3_Mutect_Vardict_ponChange.R
    chmod u+x /usr/local/annotate/annotate_CH_pd_docker3_Mutect_Vardict_ponChange.R
    
    echo "Value of mutect_vcf: '$mutect_vcf'"
    echo "Value of vardict_vcf: '$vardict_vcf'"
    echo "Value of tumor_sample_name: '$tumor_sample_name'"
    echo "Value of impact_annotation: '$impact_annotation'"
    echo "Value of intervals: '$intervals'"
    echo "Value of dockerimage_annotate_PD: '$dockerimage_annotate_PD'"
    echo "Value of project: '$project'"
    echo "Value of bqsr_spark: '$bqsr_spark'" # bqsr_spark job name
    echo "Value of mutect: '$mutect'" # mutect job name
    echo "Value of vardict: '$vardict'" # vardict job name
    echo "Value of cosmic_files: '${cosmic_files[@]}'"
    echo "Value of bolton_bick_vars: '$bolton_bick_vars'" # vardict job name
    echo "Value of mut1_bick: '$mut1_bick'" # vardict job name
    echo "Value of mut1_kelly: '$mut1_kelly'" # vardict job name

    dx-download-all-inputs --parallel

    if [ -z "$tumor_sample_name" ]
    then
        echo "\$tumor_sample_name is empty"
        tumor_sample_name=$(/usr/bin/bcftools query -l $mutect_vcf_path)
        echo $tumor_sample_name
    fi

    for cosmic_file in "${cosmic_files_path[@]}"; do 
        mv $cosmic_file .  
    done
    
    ## get <eid>_23153_0_0 from <eid>_23153_0_0.<suffix>
    eid_nameroot=$(echo $mutect_vcf_name | cut -d'.' -f1)
    target_length=$(grep -v "^@" $intervals_path | awk '{print $1,$2,$3}' | awk '$4=$3-$2+1 {sum+=$4}END {print sum}')

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
        -e $eid_nameroot \
        --bolton-bick-vars $bolton_bick_vars_path \
        --mut2-bick $mut2_bick_path \
        --mut2-kelly $mut2_kelly_path \
        --matches2 $muts2_path
    
    final_tsv=$(dx upload $eid_nameroot.final.tsv --brief)
    column_check=$(dx upload $eid_nameroot.columns.txt --brief)

    dx-jobutil-add-output final_tsv "$final_tsv" --class=file
    dx-jobutil-add-output column_check "$column_check" --class=file


}
