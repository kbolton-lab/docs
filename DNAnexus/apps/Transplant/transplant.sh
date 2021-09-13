function TRANSPLANT {
	ls /Users/brian/Bolton/UKBB/clinicaldata_ukbb_transplant.csv
	tail -n+2 /Users/brian/Bolton/UKBB/clinicaldata_ukbb_transplant.csv | awk -F, '{print $1}'

	for sample in $(tail -n+2 /Users/brian/Bolton/UKBB/clinicaldata_ukbb_transplant.csv | awk -F, '{print $1}' | head -n10 )
	
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	subfolder=${1:0:2}
	project=$2
	echo $1 $subfolder $wf_folder $project
	/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template3.sh $1 $subfolder $wf_folder $project

	echo dx run workflow-G4xzYxjJQ28FFXY59fv1ZxK7 -f $wf_folder/input_json/$subfolder/$1.json -y --priority low
	echo $1 >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log12
}