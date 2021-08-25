eid=$1
folder_prefix=$2 # i.e. 10-60
wf_folder=$3
cram=$(dx find data --class file --path UKBB_Exome_2021:"/Bulk/Exome sequences/Exome OQFE CRAM files/$folder_prefix" --json --name ${eid}_23153_0_0.cram | jq -c '.[].describe.id')
cram_index=$(dx find data --class file --path UKBB_Exome_2021:"/Bulk/Exome sequences/Exome OQFE CRAM files/$folder_prefix" --json --name ${eid}_23153_0_0.cram.crai | jq -c '.[].describe.id')

/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/cram_in2.py $eid $cram $cram_index $wf_folder

#/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/Mutect_Vardict_CH_NO_BQSR


