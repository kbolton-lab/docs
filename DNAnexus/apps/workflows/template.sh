eid=$1
file=$(dx find data --class file --path UKBB_Exome_2021:"/Bulk/Exome sequences/Exome OQFE CRAM files" --json --name ${eid}_23153_0_0.cram | jq -c '.[].describe.id')

out_folder=$3
/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/cram_in.py $eid $file $2

# /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/Mutect_Vardict_CH_NO_BQSR



