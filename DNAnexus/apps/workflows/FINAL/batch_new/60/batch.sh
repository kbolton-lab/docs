folder=60

# dx ls "$NEW_EXOME:Bulk/Exome sequences/Exome OQFE CRAM files/$folder/*.cram" | cut -d_ -f1 | sort | uniq > samples.txt
# PROJECT=$NEW_EXOME

dx build --workflow --destination $NEW_EXOME:Brian/workflows/ .
workflow=workflow-GG3zpj0JP6z34Jqz197FYKpZ


dx generate_batch_inputs -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam='(.*)_23143_0_0.cram$' -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index='(.*)_23143_0_0.cram.crai' --path "$NEW_EXOME:Bulk/Exome sequences/Exome OQFE CRAM files/$folder" -o "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch"
dx find executions -n5000 --project $NEW_EXOME --state in_progress --origin-jobs | cut -d ' ' -f2 | grep -v "^$" | cut -d- -f2 | sort > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/in_progress.txt
dx ls "project-G8YxJ2QJP6z69ZVVFVfBYyPB:/Brian/Workflow_Outputs/MutectVEP/$folder/*.vcf.gz" | cut -d_ -f1 | sort > mutect.germline
dx ls "project-G8YxJ2QJP6z69ZVVFVfBYyPB:/Brian/Workflow_Outputs/VardictIntersectMutect/$folder/*.vcf.gz" | cut -d_ -f1 | sort > vardict.germline
grep -f mutect.germline vardict.germline > both.txt
for batch_file in $(ls dx_batch.00**.tsv); do num=$(echo $batch_file | cut -d. -f2); cat both.txt samples.original.txt /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/in_progress.txt | grep -v -f /dev/stdin dx_batch.$num.tsv > dx_batch.$num.tsv.tmp && mv dx_batch.$num.tsv.tmp dx_batch.$num.tsv; done
wc -l dx_batch.00*.tsv


########################################################################################################
########################################################################################################
########################################################################################################
wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.*.tsv
## run all since not many
for batch in {0000..004}; do
   word_count=$(wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.$batch.tsv | awk '{print $1}')
   if [ $word_count -gt 1 ]; then
      dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.$batch.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 4=mem2_ssd1_v2_x16
   else
   echo no
   fi
done








# 0000
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0000.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 2=mem2_ssd1_v2_x4

# 0001
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0001.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem1_ssd1_v2_x16 --instance-type 5=mem1_ssd1_v2_x16 --instance-type 2=mem2_ssd1_v2_x4

# 0002
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0002.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}'  --instance-type 2=mem2_ssd1_v2_x4

# 0003
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0003.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem2_ssd2_x8 --instance-type 2=mem2_ssd1_v2_x4

# 0004
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0004.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'  --instance-type 2=mem2_ssd1_v2_x4 --instance-type 4=mem2_ssd1_v2_x16