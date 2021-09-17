function TRANSPLANT {
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	project=$KELLY
	# tail -n+2 /Users/brian/Bolton/UKBB/clinicaldata_ukbb_transplant.csv | awk -F, '{print $1}' | grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/total.txt

	i=0
	for sample in $(grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/total.txt /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/Transplant/eids.txt | tail -n+23); do
		subfolder=${sample:0:2}
		/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template3.sh $sample $subfolder $wf_folder $project
		dx run workflow-G500ZX0JQ28B6YFXPk32842P -f $wf_folder/input_json/$subfolder/$sample.json -y --priority low
		echo $sample >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log${subfolder}
		i=$((i+1))
		if [[ $i -eq 80 || $i -eq 160 || $i -eq 240 || $i -eq 320 ]]
		then
			echo "YES"
			sleep 1200
		fi
	done
}