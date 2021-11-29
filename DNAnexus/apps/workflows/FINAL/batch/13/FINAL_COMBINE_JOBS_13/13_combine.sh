workflow_id=workflow-G53QpP0J6XGGK9Xj6Qf9gGBv
folder=13
dx generate_batch_inputs -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam='(.*)_23153_0_0.cram$' -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index='(.*)_23153_0_0.cram.crai' --path "UKBB_Exome_2021:Bulk/Exome sequences/Exome OQFE CRAM files/$folder" -o "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch"
PROJECT=$BRIAN

for folder in 13; do
    cd /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder
    for batch in dx_batch.0000.tsv; do
        ls $batch
        tail -n +2 $batch | split -l 100 - split_${batch}_
        for file in split_${batch}_*; do
            head -n 1 $batch > ${batch}_tmp_file
            cat $file >> ${batch}_tmp_file
            mv -f ${batch}_tmp_file "$file"
            wc -l $file
            dx run $workflow_id --batch-tsv $file -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$PROJECT -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}, "stageSystemRequirements": {"stage-G4KqPz0J6XGJZGB842qJVYQK":{"executionPolicy": {"restartOn": {"AppInternalError": 1}}}}}'
            sleep 1200
        done
    done
done


