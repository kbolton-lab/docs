dx run workflow-G6Xfjj0JQ28P24VG0p41ZpF3 --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/test.hc_caller.batch -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'




workflow=workflow-G6Xfjj0JQ28P24VG0p41ZpF3
function GERMLINE_TP53 {
	i=1
	batch=$1
	workflow=$2
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.80
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.160
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.240
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.320
	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.400

	batch_total=$(wc -l $batch | awk '{print $1}')
	for sample in $(cat $batch); do
		subfolder=${sample:0:2}
		cram=$(dx find data --class file --path UKBB_Exome_2021:"/Bulk/Exome sequences/Exome OQFE CRAM files/$subfolder" --json --name ${sample}_23153_0_0.cram | jq -r '.[].describe.id')
		cram_index=$(dx find data --class file --path UKBB_Exome_2021:"/Bulk/Exome sequences/Exome OQFE CRAM files/$subfolder" --json --name ${sample}_23153_0_0.cram.crai | jq -r '.[].describe.id')
		echo $sample
		echo cram=$cram
		echo cram_index=$cram_index
		
		if [[ $i -lt 81 ]]; then
			printf "$sample	${sample}_23153_0_0.cram\t${sample}_23153_0_0.cram.crai\tproject-G3Yj1vjJ6XG579jbKyjXPGGY:$cram\tproject-G3Yj1vjJ6XG579jbKyjXPGGY:$cram_index\n" >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.80
		elif [[ $i -lt 161 ]]; then
			printf "$sample	${sample}_23153_0_0.cram	${sample}_23153_0_0.cram.crai	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram_index\n" >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.160
		elif [[ $i -lt 241 ]]; then
			printf "$sample	${sample}_23153_0_0.cram	${sample}_23153_0_0.cram.crai	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram_index\n" >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.240
		elif [[ $i -lt 321 ]]; then
			printf "$sample	${sample}_23153_0_0.cram	${sample}_23153_0_0.cram.crai	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram_index\n" >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.320
		else
			printf "$sample	${sample}_23153_0_0.cram	${sample}_23153_0_0.cram.crai	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram	project-G3Yj1vjJ6XG579jbKyjXPGGY:$cram_index\n" >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.400
		fi

		if [[ $batch_total -lt 80 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.80 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 80  ]]; 
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.80 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
			sleep 1
		elif [[ $batch_total -lt 160 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.160 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 160 ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.160 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
			sleep 1
		elif [[ $batch_total -lt 240 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.240 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 240 ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.240 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
			sleep 1
		elif [[ $batch_total -lt 320 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.320 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 320 ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.320 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
			sleep 1
		elif [[ $batch_total -lt 400 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.400 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 400 ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.400 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		fi
		i=$((i+1))
	done
}
export -f GERMLINE_TP53
GERMLINE_TP53 /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/germline.bqsr.ran.txt $workflow


workflow=workflow-G6Xj1F0JQ28J2G21FqkjgKJ8
function GERMLINE_TP53 {
	i=1
	batch=$1
	workflow=$2

	cp /Users/brian/Bolton/UKBB/batch.header /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/splits/$(basename $batch).batch.header.400

	batch_total=$(wc -l $batch | awk '{print $1}')
	for sample in $(cat $batch); do
		subfolder=${sample:0:2}
		cram=$(dx find data --class file --path UKBB_Exome_2021:"/Bulk/Exome sequences/Exome OQFE CRAM files/$subfolder" --json --name ${sample}_23153_0_0.cram | jq -r '.[].describe.id')
		cram_index=$(dx find data --class file --path UKBB_Exome_2021:"/Bulk/Exome sequences/Exome OQFE CRAM files/$subfolder" --json --name ${sample}_23153_0_0.cram.crai | jq -r '.[].describe.id')
		echo $sample
		echo cram=$cram
		echo cram_index=$cram_index
		
		
        printf "$sample	${sample}_23153_0_0.cram\t${sample}_23153_0_0.cram.crai\tproject-G3Yj1vjJ6XG579jbKyjXPGGY:$cram\tproject-G3Yj1vjJ6XG579jbKyjXPGGY:$cram_index\n" >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/splits/$(basename $batch).batch.header.400
		

		if [[ $batch_total -lt 400 && $i -eq $batch_total ]];
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/splits/$(basename $batch).batch.header.400 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
		elif [[ $i -eq 400  ]]; 
		then
			echo "YES - $i"
			echo dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/$(basename $batch).batch.header.400 -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$project -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
			sleep 1
		fi
		i=$((i+1))
	done
}
export -f GERMLINE_TP53
GERMLINE_TP53 /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/GERMLINE_TP53/splits/split_aa $workflow