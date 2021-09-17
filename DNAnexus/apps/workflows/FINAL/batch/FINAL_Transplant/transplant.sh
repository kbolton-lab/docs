function TRANSPLANT {
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	project=$KELLY


	i=0
	for sample in $(grep -v 5977335 /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/Transplant/eids.txt | sort | uniq); do
		subfolder=${sample:0:2}
		/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template3.sh $sample $subfolder $wf_folder $project
		dx run workflow-G4z4k0jJQ2888Xxv7PFkPZk2 -f $wf_folder/input_json/$subfolder/$sample.json -y --priority low --project $KELLY
		echo $sample >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log${subfolder}
		i=$((i+1))
		if [[ $i -eq 80 || $i -eq 160 || $i -eq 240 || $i -eq 320 ]]
		then
			echo "YES"
			sleep 12
		fi
	done
}
# sample=5977335
#/CH_Exome/Workflow_Outputs/FINAL/VardictJava/Transplant
#/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/Transplant
#/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Transplant