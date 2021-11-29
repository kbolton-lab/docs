

workflow=workflow-G6F0yKQJQ281QJGpKzzXQQGQ
function GERMLINE {
	i=1
	batch=$1
	workflow=$2
	project=$3
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.80
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.160
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.240
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.320
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.400

	batch_total=$(wc -l $batch | awk '{print $1}')
	for sample in $(cat $batch); do
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
			sleep 1
		elif [[ $batch_total -lt 160 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.160 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 160 ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.160 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
			sleep 1
		elif [[ $batch_total -lt 240 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.240 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 240 ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.240 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
			sleep 1
		elif [[ $batch_total -lt 320 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.320 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 320 ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/$(basename $batch).batch.header.320 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
			sleep 1
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
export -f GERMLINE 


dx ls project-G3Yj1vjJ6XG579jbKyjXPGGY:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Germline/"*.tsv" | cut -d_ -f1 | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/brian.complete.txt
for batch in aa ac ae ag ai ak am ao aq as au aw; do
    grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/brian.complete.txt /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_${batch} > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_${batch}_remaining
done
wc -l  /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_a{a,c,e,g,i,k,m,o,q,s,u,w}_remaining
for batch in aa ac ae ag ai ak am ao aq as au aw; do
    GERMLINE /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_${batch}_remaining $workflow $BRIAN
done
wc -l  /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_a{a,c,e,g,i,k,m,o,q,s,u,w}_remaining

dx run workflow-G6F0yKQJQ281QJGpKzzXQQGQ --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_as_remaining.batch.header.80  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$BRIAN -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
dx run workflow-G6F0yKQJQ281QJGpKzzXQQGQ --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_au_remaining.batch.header.80  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$BRIAN -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
dx run workflow-G6F0yKQJQ281QJGpKzzXQQGQ --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_au_remaining.batch.header.160  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$BRIAN -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
dx run workflow-G6F0yKQJQ281QJGpKzzXQQGQ --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_au_remaining.batch.header.240  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$BRIAN -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
dx run workflow-G6F0yKQJQ281QJGpKzzXQQGQ --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_au_remaining.batch.header.320  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$BRIAN -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'



dx ls project-G4qpk1jJQ285yvbXPFZKXkk8:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Germline/"*.tsv" | cut -d_ -f1 | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/kelly.complete.txt
for batch in ab ad af ah aj al an ap ar at av; do
    grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/kelly.complete.txt /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_${batch} > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_${batch}_remaining
done
for batch in ab ad af ah aj al an ap ar at av; do
    GERMLINE /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_${batch}_remaining $workflow $KELLY
done
wc -l  /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_a{b,d,f,h,j,l,n,p,r,t,v}_remaining

dx run workflow-G6F0yKQJQ281QJGpKzzXQQGQ --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_ap_remaining.batch.header.80  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$KELLY -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'





for i in 80 160; do
	dx run workflow-G6F0yKQJQ281QJGpKzzXQQGQ --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_ap_remaining.batch.header.$i  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=project-G4qpk1jJQ285yvbXPFZKXkk8 -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
	sleep 1800
done
dx run workflow-G6F0yKQJQ281QJGpKzzXQQGQ --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_ar_remaining.batch.header.80  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=project-G4qpk1jJQ285yvbXPFZKXkk8 -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
dx run workflow-G6F0yKQJQ281QJGpKzzXQQGQ --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_at_remaining.batch.header.80  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=project-G4qpk1jJQ285yvbXPFZKXkk8 -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
dx run workflow-G6F0yKQJQ281QJGpKzzXQQGQ --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE/splits/split_av_remaining.batch.header.80  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=project-G4qpk1jJQ285yvbXPFZKXkk8 -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
