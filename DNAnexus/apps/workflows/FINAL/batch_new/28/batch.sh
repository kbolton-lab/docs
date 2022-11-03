folder=28

# dx ls "$NEW_EXOME:Bulk/Exome sequences/Exome OQFE CRAM files/$folder/*.cram" | cut -d_ -f1 | sort | uniq > samples.txt
# PROJECT=$NEW_EXOME

dx build --workflow --destination $NEW_EXOME:Brian/workflows/ .
workflow=workflow-GBYJ0p0JP6zP5KVK0pVYGb0y


dx generate_batch_inputs -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam='(.*)_23143_0_0.cram$' -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index='(.*)_23143_0_0.cram.crai' --path "$NEW_EXOME:Bulk/Exome sequences/Exome OQFE CRAM files/$folder" -o "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch"


dx find executions -n2000 --project $NEW_EXOME --state in_progress --origin-jobs | cut -d ' ' -f2 | grep -v "^$" | cut -d- -f2 | sort > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/in_progress.txt
folder=28
cd /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder
dx ls "project-G8YxJ2QJP6z69ZVVFVfBYyPB:/Brian/Workflow_Outputs/MutectVEP/$folder/*.vcf.gz" | cut -d_ -f1 | sort > mutect.germline
dx ls "project-G8YxJ2QJP6z69ZVVFVfBYyPB:/Brian/Workflow_Outputs/VardictIntersectMutect/$folder/*.vcf.gz" | cut -d_ -f1 | sort > vardict.germline
grep -f mutect.germline vardict.germline > both.txt
for batch_file in $(ls dx_batch.00**.tsv); do num=$(echo $batch_file | cut -d. -f2); cat both.txt samples.original.txt /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/in_progress.txt | grep -v -f /dev/stdin dx_batch.$num.tsv > dx_batch.$num.tsv.tmp && mv dx_batch.$num.tsv.tmp dx_batch.$num.tsv; done
wc -l dx_batch.00*.tsv

# tail -n+2 cleanup/Alignment_QC.txt | cut -d, -f4 | cut -d_ -f1 | sort | uniq > test.1
# grep -v -f test.1 $folder.crams.run > test.2 
# dx ls project-G3Yj1vjJ6XG579jbKyjXPGGY:/CH_Exome/Workflow_Outputs/FINAL/BQSR/$folder/ | cut -d_ -f1 | sort | uniq > test.3
# grep -f test.2 test.3 > bqsr.txt
# for sample in $(cat bqsr.txt); do
#    bam=$(dx find data --class file --path $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/BQSR/$folder/ --name "${sample}_23153_0_0.bqsr.bam" --brief | cut -d: -f2)

#    dx run $KELLY:/test/tools/qc_exome_ukbb -iomni_vcf=file-G3zPG58J6XG2z5Qb8yGfPfP2 -ipicard_metric_accumulation_level="LIBRARY" -ibait_intervals=file-G5bgvB8JQ288QVGK2JP6BxYG -iper_base_intervals=file-G5bgvB8JQ288QVGK2JP6BxYG -itarget_intervals=file-G5bgvB8JQ288QVGK2JP6BxYG -iper_target_intervals=file-G5bgvB8JQ288QVGK2JP6BxYG -ireference=file-G3jKk3QJ6XG9KY8680qg9j8b -ireference_index=file-G3jKq0jJ6XGFqZ6Z9q4gB1yg -ibam=$bam -y --destination $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/$folder/ --name "${sample}-QC"
# done
for folder in $(dx ls $NEW_EXOME:Brian/Workflow_Outputs); do
   for dir in {10..60}; do
      dx mkdir $NEW_EXOME:Brian/Workflow_Outputs/${folder}/$dir/
   done
done

## single
name=1010590
dx run $workflow -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam=file-G9yjf08JP6zF11QPPpzJFB6P -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index=file-G9yjf4jJP6z7F1Vk6zQjV8q1 -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 2=mem2_ssd1_v2_x4 --name "$name"

########################################################################################################
########################################################################################################
########################################################################################################
wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.*.tsv
## run all since not many
for batch in {0000..0018}; do
   word_count=$(wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.$batch.tsv | awk '{print $1}')
   if [ $word_count -gt 1 ]; then
      dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.$batch.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8 
   else
   echo no
   fi
done


"mem1_ssd1_v2_x8": 0.0528,
"mem2_ssd1_v2_x8": 0.0528,
"mem3_ssd1_v2_x8": 0.0528,

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0019.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'

# GERMLINE
folder=28
workflow=workflow-GBYJ0p0JP6zP5KVK0pVYGb0y
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.germline.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 2=mem2_ssd1_v2_x4

dx ls "project-G8YxJ2QJP6z69ZVVFVfBYyPB:/Brian/Workflow_Outputs/MutectVEP/$folder/*.vcf.gz" | cut -d_ -f1 | sort > mutect.germline
dx ls "project-G8YxJ2QJP6z69ZVVFVfBYyPB:/Brian/Workflow_Outputs/VardictIntersectMutect/$folder/*.vcf.gz" | cut -d_ -f1 | sort > vardict.germline
grep -f mutect.germline vardict.germline > both.txt
dx find executions -n500 --project $NEW_EXOME --state in_progress --origin-jobs | cut -d ' ' -f2 | grep -v "^$" | cut -d- -f2 | sort > in_progress.txt
cat both.txt in_progress.txt | grep -v -f /dev/stdin /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.germline.tsv > germline.rerun.tsv
dx run $workflow --batch-tsv germline.rerun.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'

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

# 0005
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0005.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem1_ssd1_v2_x8

# 0006
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0006.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 2=mem2_ssd1_v2_x4

# 0007
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0007.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem1_ssd1_v2_x16 --instance-type 5=mem1_ssd1_v2_x16 --instance-type 2=mem2_ssd1_v2_x4

# 0008
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0008.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}'  --instance-type 2=mem2_ssd1_v2_x4
sleep 3600
# 0009
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0009.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem2_ssd2_x8 --instance-type 2=mem2_ssd1_v2_x4

# 0010
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0010.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'  --instance-type 2=mem2_ssd1_v2_x4 --instance-type 4=mem2_ssd1_v2_x16

# 0011
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0011.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem1_ssd1_v2_x8
sleep 3600
# 0012
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0012.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 2=mem2_ssd1_v2_x4

# 0013
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0013.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem1_ssd1_v2_x16 --instance-type 5=mem1_ssd1_v2_x16 --instance-type 2=mem2_ssd1_v2_x4

# 0014
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0014.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 2=mem2_ssd1_v2_x4
sleep 3600
# 0015
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0015.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem2_ssd2_x8 --instance-type 2=mem2_ssd1_v2_x4

# 0016
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0016.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 2=mem2_ssd1_v2_x4 --instance-type 4=mem2_ssd1_v2_x16

# 0017
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0017.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem1_ssd1_v2_x8
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0018.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem1_ssd1_v2_x8
#Any in progress?


cp dxworkflow.json ../$((folder + 1))/dxworkflow.json && gsed -i "s|/$folder|/$((folder + 1))|g" ../$((folder + 1))/dxworkflow.json

for i in {17..60}; do;
   cp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/batch.sh  /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$i
   gsed -i "s/folder=28/folder=$i/g" /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$i/batch.sh
   cp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dxworkflow.json  /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$i
   gsed -i "s/\"name\": \"14\"/\"name\": \"$i\"/" /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$i/dxworkflow.json
   gsed -i "s/Batch_New_14/Batch_New_${i}/" /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$i/dxworkflow.json
   gsed -i "s|/14|/$i|g" /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$i/dxworkflow.json
done

workflow=workflow-G9yjkV0JP6z293px09Xx9V7V