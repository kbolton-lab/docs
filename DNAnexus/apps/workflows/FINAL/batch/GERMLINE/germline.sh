dx login
workflow=workflow-G6B7FY8JQ28FqzfQP1pzv62f
function GERMLINE_K {
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
			sleep 1200
		fi
	done
}
export -f GERMLINE_K
GERMLINE_K /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_ar $workflow $KELLY



dx login 
workflow=workflow-G6B7FY8JQ28FqzfQP1pzv62f
function GERMLINE_B {
	i=1
	batch=$1
	workflow=$2
	project=$3
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.80
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.160
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.240
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.320
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.400
	i=81
	batch_total=$(wc -l $batch | awk '{print $1}')
	for sample in $(cat $batch | tail -n+81); do
		subfolder=${sample:0:2}
		cram=$(dx find data --class file --path UKBB_Exome_2021:"/Bulk/Exome sequences/Exome OQFE CRAM files/$subfolder" --json --name ${sample}_23153_0_0.cram | jq -r '.[].describe.id')
		cram_index=$(dx find data --class file --path UKBB_Exome_2021:"/Bulk/Exome sequences/Exome OQFE CRAM files/$subfolder" --json --name ${sample}_23153_0_0.cram.crai | jq -r '.[].describe.id')
		echo $sample
		echo cram=$cram
		echo cram_index=$cram_index
		
		if [[ $i -lt 81 ]]; then
			printf "$sample	${sample}_23153_0_0.cram\t${sample}_23153_0_0.cram.crai\tproject-G3Yj1vjJ6XG579jbKyjXPGGY:$cram\tproject-G3Yj1vjJ6XG579jbKyjXPGGY:$cram_index\n" >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.80
		elif [[ $i -lt 161 ]]; then
			printf "$sample	${sample}_23153_0_0.cram	${sample}_23153_0_0.cram.crai	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram_index\n" >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.160
		elif [[ $i -lt 241 ]]; then
			printf "$sample	${sample}_23153_0_0.cram	${sample}_23153_0_0.cram.crai	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram_index\n" >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.240
		elif [[ $i -lt 321 ]]; then
			printf "$sample	${sample}_23153_0_0.cram	${sample}_23153_0_0.cram.crai	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram_index\n" >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.320
		else
			printf "$sample	${sample}_23153_0_0.cram	${sample}_23153_0_0.cram.crai	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram_index\n" >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.400
		fi
		
		echo $sample >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log${subfolder}

		if [[ $batch_total -lt 80 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.80 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 80  ]]; 
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.80 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
			sleep 30
		elif [[ $batch_total -lt 160 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.160 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 160 ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.160 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
			sleep 30
		elif [[ $batch_total -lt 240 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.240 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 240 ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.240 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
			sleep 30
		elif [[ $batch_total -lt 320 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.320 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 320 ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.320 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
			sleep 30
		elif [[ $batch_total -lt 400 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.400 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 400 ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.400 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		fi
		i=$((i+1))
	done
}
export -f GERMLINE_B
GERMLINE_B /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_aq $workflow $BRIAN
GERMLINE_B /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_ar $workflow $KELLY

function GERMLINE_B {
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
			sleep 1200
		fi
	done
}
export -f GERMLINE_B
GERMLINE_B /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_aq $workflow $BRIAN


# single sample
dx login 
workflow=workflow-G6B7FY8JQ28FqzfQP1pzv62f
function GERMLINE_S {
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
export -f GERMLINE_S
GERMLINE_S $workflow 1281328 $KELLY





for sample in $(cat /Users/brian/Bolton/UKBB/results/lung/update/brian.redo.txt | tail -n+2); do
	GERMLINE_S $workflow $sample $BRIAN
done

for sample in $(cat /Users/brian/Bolton/UKBB/results/lung/update/kelly.redo.txt | tail -n+3); do
	GERMLINE_S $workflow $sample $KELLY
done

for sample in $(cat kelly.redo.txt | head); do
	#wc -l ${sample}_23153_0_0.final.tsv #> 
	wc -l ${sample}_23153_0_0.pon.pass.tsv
done

for sample in $(cat brian.redo.txt); do
	awk -F'\t' '$202 < 1.260958e-09 {print}' OFS='\t' ${sample}_23153_0_0.final.tsv > ${sample}_23153_0_0.pon.pass.tsv
done