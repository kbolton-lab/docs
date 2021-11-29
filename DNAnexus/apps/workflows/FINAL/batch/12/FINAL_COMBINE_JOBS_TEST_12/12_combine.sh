project=$KELLY
workflow_id=workflow-G53GjJ0JQ2811zB60BPf6pB0

for folder in 12; do
    cd /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder
    dx find analyses --created-after 2021-09-19 -n5000 --state failed --name "FINAL_COMBINE_JOBS_$folder-*" --project=$PROJECT | grep failed | cut -d' ' -f2 | cut -d- -f2 > FINAL_COMBINE_JOBS_$folder/remain_failed.txt 
    dx find analyses --created-after 2021-09-19 -n5000 --state terminated  --name "FINAL_COMBINE_JOBS_$folder-*" --project=$PROJECT  | grep terminated | cut -d' ' -f2 | cut -d- -f2  > FINAL_COMBINE_JOBS_$folder/remain_terminated.txt
    for batch in dx_batch.000{3,4,5}; do
        (head -n1 $batch.tsv; grep -f FINAL_COMBINE_JOBS_12/remain_failed.txt $batch.tsv) > $batch.transplant.match.tsv
        # grep -v -f FINAL_COMBINE_JOBS_12/remaining.transplant.match.txt $batch.tsv > $batch.latest.tsv
        dx run $workflow_id --batch-tsv $batch.transplant.match.tsv -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$KELLY -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}, "stageSystemRequirements": {"stage-G4KqPz0J6XGJZGB842qJVYQK":{"executionPolicy": {"restartOn": {"AppInternalError": 1}}}}}' 
        sleep 1
    done
done

wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/dx_batch.000{3,4,5}.trans*.tsv
323 + 3 = 326
wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/dx_batch.000{3,4,5}.latest.tsv
1177+3=1180
stage-G4KqPz0J6XGJZGB842qJVYQK

dx find analyses --created-after 2021-09-19 --state "done" --name "FINAL_COMBINE_JOBS_12-*" --project $KELLY