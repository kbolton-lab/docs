dx generate_batch_inputs -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam='(.*)_23153_0_0.cram$' -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index='(.*)_23153_0_0.cram.crai' --path "UKBB_Exome_2021:Bulk/Exome sequences/Exome OQFE CRAM files/11"

cp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/test/FINAL/dxworkflow.json FINAL_11
gsed -i 's|/10|/11|g' FINAL_11/dxworkflow.json
dx build FINAL_11 --workflow --destination /test/workflows/ --keep-open
{"id": "workflow-G4gzV7QJ6XG2Q4yX4K4g4FPQ"}

echo '
function W11_FINAL {
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template2.sh $1 $2 $wf_folder

	dx run workflow-G4gzV7QJ6XG2Q4yX4K4g4FPQ -f $wf_folder/input_json/$2/$1.json -y
	echo $1 >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log11
}' >> 11.sh

grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log11 dx_batch.0000.tsv > dx_batch.0000.tsv.tmp && mv dx_batch.0000.tsv.tmp dx_batch.0000.tsv
dx run workflow-G4gzV7QJ6XG2Q4yX4K4g4FPQ --batch-tsv dx_batch.0000.tsv --name FINAL_11 -y --priority low
