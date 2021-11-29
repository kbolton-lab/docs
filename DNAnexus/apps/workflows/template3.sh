eid=$1
folder_prefix=$2 # i.e. 10-60
wf_folder=$3
project=$4
cram=$(dx find data --class file --path UKBB_Exome_2021:"/Bulk/Exome sequences/Exome OQFE CRAM files/$folder_prefix" --json --name ${eid}_23153_0_0.cram | jq -c '.[].describe.id')
cram_index=$(dx find data --class file --path UKBB_Exome_2021:"/Bulk/Exome sequences/Exome OQFE CRAM files/$folder_prefix" --json --name ${eid}_23153_0_0.cram.crai | jq -c '.[].describe.id')
echo cram=$cram
echo cram_index=$cram_index

echo $eid $cram $cram_index $wf_folder $folder_prefix $project 
/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/cram_in3.py $eid $cram $cram_index $wf_folder $folder_prefix $project



