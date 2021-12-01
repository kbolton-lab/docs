folder=12
folder_m1=$(($folder - 1 ))
cp -r /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder_m1 /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder
mv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/$folder_m1.sh mv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/$folder.sh

mv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/FINAL_${folder_m1} /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/FINAL_${folder}
gsed -i "s/FINAL_${folder_m1}/FINAL_${folder}/g" /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/FINAL_${folder}/dxworkflow.json

gsed -i "s|/${folder_m1}|/${folder}|g" /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/FINAL_${folder}/dxworkflow.json
workflow_id=$(dx build /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/FINAL_${folder} --workflow --destination /test/workflows/ --keep-open | jq -r .id )

#"id": "workflow-G4qkB9QJ6XGKVJJ1Bxv6Gz2y"
mv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/${folder_m1}.sh /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/${folder}.sh
gsed -i "s/W${folder_m1}_FINAL/W${folder}_FINAL/g" /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/${folder}.sh
gsed -i "s|logs/log${folder_m1}|logs/log${folder}|g" /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/${folder}.sh
gsed -i "s/workflow-[a-zA-Z0-9]\+/$workflow_id/g" /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/${folder}.sh
touch /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log12

rm /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch*
dx generate_batch_inputs -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam='(.*)_23153_0_0.cram$' -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index='(.*)_23153_0_0.cram.crai' --path "UKBB_Exome_2021:Bulk/Exome sequences/Exome OQFE CRAM files/$folder" -o "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch"
source /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/${folder}.sh
mkdir /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/input_json/$folder
1=$(sed -n '2p' /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0000.tsv | awk '{print $1}')
2=$folder
W${folder}_FINAL $1 $2

grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/logs/log${folder} /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0000.tsv > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0000.tsv.tmp && mv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0000.tsv.tmp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0000.tsv

# workflow_id=workflow-G4qkQJQJ6XG0P2Bz4Qj3zxk2
# dx run $workflow_id --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0000.tsv -y --priority low

<<<<<<< HEAD
<<<<<<< HEAD
workflow_id=workflow-G4x6JB0JQ281gy5Q4K814GXF
=======
workflow_id=workflow-G4x6JB0JQ281gy5Q4K814GXF # Kelly
>>>>>>> 65d148f5d1ffac4ed8440f892e4a8b651eb53b01
=======
workflow_id=workflow-G4x6JB0JQ281gy5Q4K814GXF # Kelly
>>>>>>> 65d148f5d1ffac4ed8440f892e4a8b651eb53b01
for folder in 12; do
    cd /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder
    for batch in dx_batch.0001.tsv; do
        ls $batch
        tail -n +2 $batch | split -l 100 - split_${batch}_
        for file in split_${batch}_*; do
            head -n 1 $batch > ${batch}_tmp_file
<<<<<<< HEAD
            # cat "$file" >> ${batch}_tmp_file
            grep -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/dx_batch_1_latest.0001.tsv "$file" >> ${batch}_tmp_file
            mv -f ${batch}_tmp_file "$file"
<<<<<<< HEAD
            wc -l $file
            dx run $workflow_id --batch-tsv $file -istage-G4fq8KjJ6XG6pBYqK1VXFgZV.project=$KELLY -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}'
            sleep 120
=======
            echo dx run $workflow_id --batch-tsv $file -istage-G4fq8KjJ6XG6pBYqK1VXFgZV.project=$KELLY -y --priority low
            
>>>>>>> 65d148f5d1ffac4ed8440f892e4a8b651eb53b01
=======
            cat "$file" >> ${batch}_tmp_file
            mv -f ${batch}_tmp_file "$file"
            echo dx run $workflow_id --batch-tsv $file -istage-G4fq8KjJ6XG6pBYqK1VXFgZV.project=$KELLY -y --priority low
            
>>>>>>> 65d148f5d1ffac4ed8440f892e4a8b651eb53b01
        done
    done
done

<<<<<<< HEAD
<<<<<<< HEAD
workflow_id=workflow-G4x8k0QJ6XG88z719jyVVX77
=======
workflow_id=workflow-G4x8k0QJ6XG88z719jyVVX77 # Brian
>>>>>>> 65d148f5d1ffac4ed8440f892e4a8b651eb53b01
=======
workflow_id=workflow-G4x8k0QJ6XG88z719jyVVX77 # Brian
>>>>>>> 65d148f5d1ffac4ed8440f892e4a8b651eb53b01
for folder in 12; do
    cd /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder
    for batch in dx_batch.0002.tsv; do
        ls $batch
        tail -n +2 $batch | split -l 100 - split_${batch}_
        for file in split_${batch}_*; do
            head -n 1 $batch > ${batch}_tmp_file
<<<<<<< HEAD
<<<<<<< HEAD
            grep -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/dx_batch_1_latest.0002.tsv "$file" >> ${batch}_tmp_file
            mv -f ${batch}_tmp_file "$file"
            wc -l $file
            dx run $workflow_id --batch-tsv $file -istage-G4fq8KjJ6XG6pBYqK1VXFgZV.project=$BRIAN -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}, "stageSystemRequirements": {"stage-G4KqPz0J6XGJZGB842qJVYQK":{"executionPolicy": {"restartOn": {"AppInternalError": 1}}}}}'
            sleep 30
=======
=======
>>>>>>> 65d148f5d1ffac4ed8440f892e4a8b651eb53b01
            cat "$file" >> ${batch}_tmp_file
            mv -f ${batch}_tmp_file "$file"
            dx run $workflow_id --batch-tsv $file -istage-G4fq8KjJ6XG6pBYqK1VXFgZV.project=$BRIAN -y --priority low
            
<<<<<<< HEAD
>>>>>>> 65d148f5d1ffac4ed8440f892e4a8b651eb53b01
=======
>>>>>>> 65d148f5d1ffac4ed8440f892e4a8b651eb53b01
        done
    done
done