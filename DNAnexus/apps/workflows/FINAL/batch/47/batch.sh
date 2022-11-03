# Brian's Project: $BRIAN
# ODDS

# Kelly's Project: $KELLY
# EVENS

########################################################################################################
########################################################################################################

folder=47
PROJECT=$BRIAN
workflow=workflow-G7v1Fq0J6XGBGp7J3GP32fYV
dx generate_batch_inputs -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam='(.*)_23153_0_0.cram$' -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index='(.*)_23153_0_0.cram.crai' --path "UKBB_Exome_2021:Bulk/Exome sequences/Exome OQFE CRAM files/$folder" -o "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch"

tail -n+2 cleanup/Alignment_QC.txt | cut -d, -f4 | cut -d_ -f1 | sort | uniq > test.1
grep -v -f test.1 $folder.crams.run > test.2 
dx ls project-G3Yj1vjJ6XG579jbKyjXPGGY:/CH_Exome/Workflow_Outputs/FINAL/BQSR/$folder/ | cut -d_ -f1 | sort | uniq > test.3
grep -f test.2 test.3 > bqsr.txt
for sample in $(cat bqsr.txt); do
   bam=$(dx find data --class file --path $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/BQSR/$folder/ --name "${sample}_23153_0_0.bqsr.bam" --brief | cut -d: -f2)

   dx run $KELLY:/test/tools/qc_exome_ukbb -iomni_vcf=file-G3zPG58J6XG2z5Qb8yGfPfP2 -ipicard_metric_accumulation_level="LIBRARY" -ibait_intervals=file-G5bgvB8JQ288QVGK2JP6BxYG -iper_base_intervals=file-G5bgvB8JQ288QVGK2JP6BxYG -itarget_intervals=file-G5bgvB8JQ288QVGK2JP6BxYG -iper_target_intervals=file-G5bgvB8JQ288QVGK2JP6BxYG -ireference=file-G3jKk3QJ6XG9KY8680qg9j8b -ireference_index=file-G3jKq0jJ6XGFqZ6Z9q4gB1yg -ibam=$bam -y --destination $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/$folder/ --name "${sample}-QC"
done


wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.*.*.tsv
########################################################################################################
########################################################################################################
########################################################################################################
## run all since not many
for batch in {0000..0008}; do
   for file in a b; do
      word_count=$(wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.$batch.$file.tsv | awk '{print $1}')
      if [ $word_count -gt 1 ]; then
         dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.$batch.$file.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem1_ssd1_v2_x8
      else
      echo no
      fi
   done
done


"mem1_ssd1_v2_x8": 0.0528,
"mem2_ssd1_v2_x8": 0.0528,
"mem3_ssd1_v2_x8": 0.0528,

# 0000
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0000.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 2=mem2_ssd1_v2_x4

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0000.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 5=mem2_ssd2_x8 

# 0001
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0001.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem1_ssd1_v2_x16 --instance-type 5=mem1_ssd1_v2_x16 --instance-type 2=mem2_ssd1_v2_x4
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0001.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  

# 0002
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0002.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}'  --instance-type 2=mem2_ssd1_v2_x4
sleep 2400
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0002.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem2_ssd1_v2_x16


# 0003

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0003.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem2_ssd2_x8 --instance-type 2=mem2_ssd1_v2_x4

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0003.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 2=mem2_ssd1_v2_x16

# 0004

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0004.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'  --instance-type 2=mem2_ssd1_v2_x4 --instance-type 4=mem2_ssd1_v2_x16
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0004.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8

# 0005

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0005.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem1_ssd1_v2_x8
sleep 1800
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0005.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8

# 0006
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0006.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'  --instance-type 5=mem2_ssd2_x8 --instance-type 2=mem2_ssd1_v2_x4

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0006.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 5=mem1_ssd1_v2_x8


# 0007

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0007.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 2=mem3_ssd1_v2_x4

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0007.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' 

# 0008
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0008.a.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'

#Any in progress?




cp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/47/batch.sh /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/49/batch.sh
cp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/47/dxworkflow.json /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/49/dxworkflow.json
