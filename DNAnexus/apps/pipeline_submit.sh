# subfolders=({10..60})

subfolder=10
samples=($(dx ls "/Bulk/Exome sequences/Exome OQFE CRAM files/$subfolder/*.cram" | cut -d'_' -f1 | sort | head -n5 | ggrep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL_WITH_BQSR2/log))

for sample in ${samples[@]}; do
    W10_EID $sample $subfolder
done


## app submit
samples=($(dx ls "/CH_Exome/Inputs/youngest/30/*.bam" | cut -d'_' -f1 | sort))

for sample in ${samples[@]:0:1}; do
    BQSR_SPARK $sample
done
