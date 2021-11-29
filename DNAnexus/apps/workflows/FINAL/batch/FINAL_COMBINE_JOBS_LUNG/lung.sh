workflow=workflow-G5qX1pjJQ28KX7qK2PXF1gpg
function LUNG_K {
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	i=0
	batch=$1
	workflow=$2
	project=$3
	for sample in $(cat $batch); do
		subfolder=${sample:0:2}
		/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template3.sh $sample $subfolder $wf_folder $project
		date
		dx run $workflow -f $wf_folder/input_json/$subfolder/$sample.json -y --priority low --project $project --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		echo $sample >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log${subfolder}
		i=$((i+1))
		if [[ $i -eq 80 || $i -eq 160 || $i -eq 240 || $i -eq 320 ]]
		then
			echo "YES"
			sleep 600
		fi
	done
}
export -f LUNG_K
LUNG_K /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/FINAL_COMBINE_JOBS_LUNG/kelly.old.lung $workflow $KELLY



dx login 

workflow=workflow-G5qX1pjJQ28KX7qK2PXF1gpg
function LUNG_B {
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	i=0
	batch=$1
	workflow=$2
	project=$3
	for sample in $(cat $batch | tail -n+2); do
		subfolder=${sample:0:2}
		/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template3.sh $sample $subfolder $wf_folder $project
		date
		dx run $workflow -f $wf_folder/input_json/$subfolder/$sample.json -y --priority low --project $project --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		
		echo $sample >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log${subfolder}
		i=$((i+1))
		if [[ $i -eq 80 || $i -eq 160 || $i -eq 240 || $i -eq 320 ]]
		then
			echo "YES"
			sleep 1200
		fi
	done
}
export -f LUNG_B
LUNG_B /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/FINAL_COMBINE_JOBS_LUNG/brian.old.lung $workflow $BRIAN


# single sample
dx login 
workflow=workflow-G5qXP6QJQ28GYkjPG9pYb8jp
function LUNG_S {
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	i=0
	workflow=$1
	sample=$2
	project=$3
	subfolder=${sample:0:2}
	/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template3.sh $sample $subfolder $wf_folder $project
	date
	dx run $workflow -f $wf_folder/input_json/$subfolder/$sample.json -y --priority low --project $project --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 3}}'
	
	echo $sample >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log${subfolder}
	i=$((i+1))
}
export -f LUNG_S
LUNG_S $workflow 1422233 $BRIAN

for sample in $(cat /Users/brian/Bolton/UKBB/results/lung/update/brian.redo.txt | tail -n+2); do
	LUNG_S $workflow $sample $BRIAN
done

for sample in $(cat /Users/brian/Bolton/UKBB/results/lung/update/kelly.redo.txt | tail -n+3); do
	LUNG_S $workflow $sample $KELLY
done

for sample in $(cat kelly.redo.txt | head); do
	#wc -l ${sample}_23153_0_0.final.tsv #> 
	wc -l ${sample}_23153_0_0.pon.pass.tsv
done

for sample in $(cat brian.redo.txt); do
	awk -F'\t' '$202 < 1.260958e-09 {print}' OFS='\t' ${sample}_23153_0_0.final.tsv > ${sample}_23153_0_0.pon.pass.tsv
done