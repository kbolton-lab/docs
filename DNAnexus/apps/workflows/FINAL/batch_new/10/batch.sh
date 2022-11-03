folder=10
#note 1067325 not in new exomes

dx ls "$NEW_EXOME:Bulk/Exome sequences/Exome OQFE CRAM files/$folder/*.cram" | cut -d_ -f1 | sort | uniq > samples.txt
PROJECT=$NEW_EXOME

dx build --workflow --destination $NEW_EXOME:Brian/workflows/ .
workflow=workflow-G9yk6g8JP6zK5X16KVy4VYVZ




dx find executions -n1000 --project $NEW_EXOME --state in_progress --origin-jobs | cut -d ' ' -f2 | grep -v "^$" | cut -d- -f2 | sort > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/in_progress.txt
folder=10
cd /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder
dx ls "project-G8YxJ2QJP6z69ZVVFVfBYyPB:/Brian/Workflow_Outputs/MutectVEP/$folder/*.vcf.gz" | cut -d_ -f1 | sort > mutect.germline
dx ls "project-G8YxJ2QJP6z69ZVVFVfBYyPB:/Brian/Workflow_Outputs/VardictIntersectMutect/$folder/*.vcf.gz" | cut -d_ -f1 | sort > vardict.germline
grep -f mutect.germline vardict.germline > both.txt
for batch_file in $(ls dx_batch.00*.tsv); do num=$(echo $batch_file | cut -d. -f2); cat both.txt samples.original.txt /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/in_progress.txt | grep -v -f /dev/stdin dx_batch.$num.tsv > dx_batch.$num.tsv.tmp && mv dx_batch.$num.tsv.tmp dx_batch.$num.tsv; done
wc -l dx_batch.00*.tsv


wc -l /Users/brian/Downloads/carrier250k_additional.txt
head /Users/brian/Downloads/carrier250k_additional.txt
t
for folder in {12..60}; do
   pfx=$(echo $folder | cut -c1,2)
   echo $pfx
   cd /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder
   grep ^$pfx /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/germline2.ids.txt > germline2.ids.$folder.txt
   wc -l dx_batch.00**.tsv | tail -n1
   wc -l germline2.ids.$folder.txt
   (head -n1 dx_batch.0000.tsv; grep -f germline2.ids.$folder.txt dx_batch.00**.tsv | cut -d: -f 2-) > dx_batch.germline2.tsv
   for batch_file in $(ls dx_batch.00*.tsv); do 
      num=$(echo $batch_file | cut -d. -f2); 
      grep -v -f germline2.ids.$folder.txt dx_batch.$num.tsv > dx_batch.$num.tsv.tmp && mv dx_batch.$num.tsv.tmp dx_batch.$num.tsv;
   done
   wc -l dx_batch.00*.tsv | tail -n1
   echo
done

# for folder in {11..16}; do
for folder in {51..60}
for folder in {42..50}
for folder in {34..41}
for folder in {25..33}
sleep 7200; sleep 7200; for folder in {51..60}; do
   wf=$(grep "^$folder" /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/folder_wfs.txt | awk '{print $2}')
   workflow=$wf
   printf "$folder\t$wf\n"
   # wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/{11..16}/dx_batch.germline2.tsv
   dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.germline2.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem1_ssd1_v2_x8
done
batch_new % cat germline*ids.txt | wc -l
36903
# not germline
dx ls "project-G8YxJ2QJP6z69ZVVFVfBYyPB:/Brian/Workflow_Outputs/BQSR/10/*.bam" | cut -d_ -f1 | sort | uniq > bqsr.txt
cat germline.ids.10.txt germline2.ids.10.txt > germline.total
grep -v -f germline.total bqsr.txt
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

name=1000128
dx run $workflow -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam=file-FyQbfk0JkF6K91ZxJ001fKY4 -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index=file-FyQbfk8JkF62jVVq8ZXGky8Q -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 2=mem2_ssd1_v2_x4 --name "$name"

########################################################################################################
########################################################################################################
########################################################################################################
wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.00*.tsv
folder=10
workflow=workflow-G9yk6g8JP6zK5X16KVy4VYVZ
## run all since not many
for batch in {0000..00018}; do
   word_count=$(wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.$batch.tsv | awk '{print $1}')
   if [ $word_count -gt 1 ]; then
      dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.$batch.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem1_ssd1_v2_x8
   else
   echo no
   fi
done


"mem1_ssd1_v2_x8": 0.0528,
"mem2_ssd1_v2_x8": 0.0528,
"mem3_ssd1_v2_x8": 0.0528,

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0019.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'

# GERMLINE
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.germline.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
dx ls "project-G8YxJ2QJP6z69ZVVFVfBYyPB:/Brian/Workflow_Outputs/MutectVEP/10/*.vcf.gz" | cut -d_ -f1 | sort > mutect.germline
dx ls "project-G8YxJ2QJP6z69ZVVFVfBYyPB:/Brian/Workflow_Outputs/VardictIntersectMutect/10/*.vcf.gz" | cut -d_ -f1 | sort > vardict.germline
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

# 0009
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0009.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem2_ssd2_x8 --instance-type 2=mem2_ssd1_v2_x4

# 0010
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0010.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'  --instance-type 2=mem2_ssd1_v2_x4 --instance-type 4=mem2_ssd1_v2_x16

# 0011
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0011.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem1_ssd1_v2_x8

# 0012
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0012.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 2=mem2_ssd1_v2_x4

# 0013
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0013.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem1_ssd1_v2_x16 --instance-type 5=mem1_ssd1_v2_x16 --instance-type 2=mem2_ssd1_v2_x4

# 0014
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0014.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 2=mem2_ssd1_v2_x4

# 0015
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0015.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem2_ssd2_x8 --instance-type 2=mem2_ssd1_v2_x4

# 0016
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0016.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 2=mem2_ssd1_v2_x4 --instance-type 4=mem2_ssd1_v2_x16

# 0017
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0017.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem1_ssd1_v2_x8

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_batch.0018.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem1_ssd1_v2_x8

#Any in progress?


cp dxworkflow.json ../$((folder + 1))/dxworkflow.json && gsed -i "s|/$folder|/$((folder + 1))|g" ../$((folder + 1))/dxworkflow.json

cp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/batch.sh /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/11/batch.sh && sed ""
dx build --workflow --destination $NEW_EXOME:Brian/workflows/ .

workflow=workflow-G9yjkV0JP6z293px09Xx9V7V