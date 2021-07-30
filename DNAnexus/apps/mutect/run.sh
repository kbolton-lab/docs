dx run /test/tools/mutect \
    -ireference="/CH_Exome/Inputs/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
    -ireference_index="/CH_Exome/Inputs/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai" \
    -ireference_dict="/CH_Exome/Inputs/GRCh38_full_analysis_set_plus_decoy_hla.dict" \
    -itumor_bam="/CH_Exome/BAMS/1010590_23153_0_0_chr22.bam" \
    -itumor_bai="/CH_Exome/BAMS/1010590_23153_0_0_chr22.bam.bai" \
    -inormal_bam="/CH_Exome/BAMS/1006108_23153_0_0_chr22.bam" \
    -inormal_bai="/CH_Exome/BAMS/1006108_23153_0_0_chr22.bam.bai" \
    -iinterval_list="/CH_Exome/Inputs/xgen_plus_spikein.GRCh38_chr22.interval_list" \
    --destination /CH_Exome/Mutect2/