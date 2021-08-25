# subfolders=({10..60})

subfolder=10
samples=($(dx ls "/Bulk/Exome sequences/Exome OQFE CRAM files/$subfolder/*.cram" | cut -d'_' -f1 | sort | head -n5 | ggrep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL_WITH_BQSR2/log))

for sample in ${samples[@]}; do
    W10_EID $sample $subfolder
done


