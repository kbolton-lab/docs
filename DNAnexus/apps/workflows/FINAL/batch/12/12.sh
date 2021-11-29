## old FINAL_12 mutect scatter
function __W12_FINAL {
	# array2=($(dx ls "/Bulk/Exome sequences/Exome OQFE CRAM files/10/" | grep -v "crai" | cut -d'_' -f1 | sort | head -n100 | grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/log))
	# for eid in ${array2[@]}; do W10_FINAL $eid 10; done
	# BQSR2
	# 1=eid
	# 2=folder_prefix=$3 # i.e. 10-60
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	2=${1:0:2}
	/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template2.sh $1 $2 $wf_folder

	# mutect _scatter
	dx run workflow-G4qkQJQJ6XG0P2Bz4Qj3zxk2 -f $wf_folder/input_json/$2/$1.json -y --priority low --project $BRIAN --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}'
	# mutect_single
	# dx run workflow-G4x6JB0JQ281gy5Q4K814GXF -f $wf_folder/input_json/$2/$1.json -y --priority low
	echo $1 >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log12
}
workflow_id=workflow-G4qkQJQJ6XG0P2Bz4Qj3zxk2
i=1
for folder in 12; do
    cd /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder
    for batch in dx_batch.0000.tsv; do
        #ls $batch
        tail -n +2 $batch | split -l 100 - split_${batch}_
        for file in split_${batch}_*; do
            head -n 1 $batch > ${batch}_tmp_file
			cat "$file" >> ${batch}_tmp_file
            mv -f ${batch}_tmp_file "$file"
            #wc -l $file
			if [[ $i -eq 1 ]]; then
				echo $i
				echo dx run $workflow_id --batch-tsv $file -istage-G4fq8KjJ6XG6pBYqK1VXFgZV.project=$BRIAN -istage-G4fq8KjJ6XG6pBYqK1VXFgZV.mutect="mutect" -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}, "stageSystemRequirements": {"stage-G4KqPz0J6XGJZGB842qJVYQK":{"executionPolicy": {"restartOn": {"AppInternalError": 1}}}}}'
			fi
            # echo dx run $workflow_id --batch-tsv $file -istage-G4fq8KjJ6XG6pBYqK1VXFgZV.project=$BRIAN -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}, "stageSystemRequirements": {"stage-G4KqPz0J6XGJZGB842qJVYQK":{"executionPolicy": {"restartOn": {"AppInternalError": 1}}}}}'
            sleep 1
			i=$(($i+1))
        done
    done
done
dx run workflow-G4qkQJQJ6XG0P2Bz4Qj3zxk2 --batch-tsv split_dx_batch.0000.tsv_aa -istage-G4fq8KjJ6XG6pBYqK1VXFgZV.project=project-G3Yj1vjJ6XG579jbKyjXPGGY -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}, "stageSystemRequirements": {"stage-G4KqPz0J6XGJZGB842qJVYQK":{"executionPolicy": {"restartOn": {"AppInternalError": 1}}}}}'




# dx run workflow-G4v1kZjJQ286B56406bBFxvQ -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/input_json/12/1208696.json -y --priority low

# kelly
# dx run workflow-G4x6JB0JQ281gy5Q4K814GXF -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/input_json/12/1219875.json -y --priority low --project project-G4qpk1jJQ285yvbXPFZKXkk8 --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}'
function W12_FINAL_K {
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	# tail -n+2 /Users/brian/Bolton/UKBB/clinicaldata_ukbb_transplant.csv | awk -F, '{print $1}' | grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/total.txt
	sample=$1
	subfolder=${sample:0:2}
	project=$2
	/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template2.sh $sample $subfolder $wf_folder $project
	echo dx run workflow-G4x6JB0JQ281gy5Q4K814GXF -f $wf_folder/input_json/$subfolder/$sample.json -y --priority low --project $KELLY --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}'
	echo $sample >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log${subfolder}
}
# dx update stage workflow-G4x6JB0JQ281gy5Q4K814GXF "mutect_vep" --executable applet-G4fpZk8J6XGGQpPQJyQXfkpQ
# dx update stage workflow-G4x6JB0JQ281gy5Q4K814GXF "final_annotation" --executable applet-G52PBG0JQ289QY9fBy9B1ZxX
# dx update stage workflow-G4x6JB0JQ281gy5Q4K814GXF "final_annotation_and_declutter" -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/FINAL_Transplant/input.json

# brian 
function W12_FINAL_B {
	wf_folder="/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL"
	# tail -n+2 /Users/brian/Bolton/UKBB/clinicaldata_ukbb_transplant.csv | awk -F, '{print $1}' | grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/total.txt
	sample=$1
	subfolder=${sample:0:2}
	project=$2
	/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/template2.sh $sample $subfolder $wf_folder $project
	
	echo dx run workflow-G4x8k0QJ6XG88z719jyVVX77 -f $wf_folder/input_json/$subfolder/$sample.json -y --priority low --project $BRIAN --extra-args '{"stageSystemRequirements": {"stage-11":{"executionPolicy": {"restartOn": {"AppInternalError": 1, "UnresponsiveWorker": 2}}}}}'
	echo $sample >> /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log${subfolder}
}

dx run workflow-G4x8k0QJ6XG88z719jyVVX77 -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/input_json/12/1233463.json -y --priority low --project project-G3Yj1vjJ6XG579jbKyjXPGGY --extra-args '{"stageSystemRequirements": {"stage-G4KqPz0J6XGJZGB842qJVYQK":{"executionPolicy": {"restartOn": {"AppInternalError": 1, "UnresponsiveWorker": 2}}}}}'

dx run workflow-G4x8k0QJ6XG88z719jyVVX77 -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/input_json/12/1237094.json -y --priority low --project project-G3Yj1vjJ6XG579jbKyjXPGGY --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}, "stageSystemRequirements": {"stage-G4KqPz0J6XGJZGB842qJVYQK":{"executionPolicy": {"restartOn": {"AppInternalError": 1}}}}}'