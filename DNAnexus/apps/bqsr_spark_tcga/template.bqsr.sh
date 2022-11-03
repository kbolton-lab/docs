eid=$1
app_folder=$2
bam=$(dx find data --class file --path UKBB_Exome_2021:"/CH_Exome/Inputs/youngest/30" --json --name ${eid}_23153_0_0.bam | jq -c '.[].describe.id')
bam_index=$(dx find data --class file --path UKBB_Exome_2021:"/CH_Exome/Inputs/youngest/30" --json --name ${eid}_23153_0_0.bam.bai | jq -c '.[].describe.id')

/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/bqsr_spark/bqsr_in.py $eid $bam $bam_index $app_folder

#/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/Mutect_Vardict_CH_NO_BQSR



