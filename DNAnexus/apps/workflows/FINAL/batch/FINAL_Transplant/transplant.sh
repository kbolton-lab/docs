function TRANSPLANT {
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	project=$KELLY
	i=0
	for sample in $(cat /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/test/failed.transplant | sort | uniq); do
		subfolder=${sample:0:2}
		/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template3.sh $sample $subfolder $wf_folder $project
		# dx run workflow-G4z4k0jJQ2888Xxv7PFkPZk2 -f $wf_folder/input_json/$subfolder/$sample.json -y --priority low --project $KELLY 
		echo dx run workflow-G4z4k0jJQ2888Xxv7PFkPZk2 -f $wf_folder/input_json/$subfolder/$sample.json -y --priority low --project $KELLY --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}'
		echo $sample >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log${subfolder}
		i=$((i+1))
		if [[ $i -eq 80 || $i -eq 160 || $i -eq 240 || $i -eq 320 ]]
		then
			echo "YES"
			sleep 12
		fi
	done
}

function TRANSPLANT {
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	project=$KELLY
	i=0
	for sample in $(cat /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/test/failed.transplant | sort | uniq | tail -n+2); do
		subfolder=${sample:0:2}
		/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template3.sh $sample $subfolder $wf_folder $project
		# dx run workflow-G4z4k0jJQ2888Xxv7PFkPZk2 -f $wf_folder/input_json/$subfolder/$sample.json -y --priority low --project $KELLY 
		dx run workflow-G4z4k0jJQ2888Xxv7PFkPZk2 -f $wf_folder/input_json/$subfolder/$sample.json -y --priority low --project $KELLY --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}'
		echo $sample >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log${subfolder}
		i=$((i+1))
		if [[ $i -eq 80 || $i -eq 160 || $i -eq 240 || $i -eq 320 ]]
		then
			echo "YES"
			sleep 12
		fi
	done
}


stage-G4KqPz0J6XGJZGB842qJVYQK.mutect_vcf
# dx update stage workflow-G4z4k0jJQ2888Xxv7PFkPZk2 "final_annotation_and_declutter" --executable applet-G52Jg78JQ287Y96V3FVXP46y
# dx update stage workflow-G4z4k0jJQ2888Xxv7PFkPZk2 "final_annotation_and_declutter" -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/FINAL_Transplant/input.json

stage-G4xz6gQJQ28JqVq0FZzPq571",
            "outputField": "annotated_pon2_vcf"

# sample=5977335
#/CH_Exome/Workflow_Outputs/FINAL/VardictJava/Transplant
#/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/Transplant
#/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Transplant

# INFO:dxpy:Deleting applet(s) applet-G51q928JQ280fVjfJjp7V9yJ
# {"id": "applet-G52Jg78JQ287Y96V3FVXP46y"}