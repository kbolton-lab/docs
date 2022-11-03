# Brian's Project: $Brian
12
13
# ODDS
15


# Kelly's Project: $KELLY
10
11
(12)
# EVENS
20 ****************************************************
22
24

folder=48
PROJECT=$KELLY
workflow=workflow-G7x2920JQ283GxJV4910qBBv
dx generate_batch_inputs -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam='(.*)_23153_0_0.cram$' -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index='(.*)_23153_0_0.cram.crai' --path "UKBB_Exome_2021:Bulk/Exome sequences/Exome OQFE CRAM files/$folder" -o "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch"

# "mem1_ssd1_v2_x8": 0.0528,
# "mem2_ssd1_v2_x8": 0.0528,
# "mem3_ssd1_v2_x8": 0.0528,


tail -n+2 cleanup/Alignment_QC.txt | cut -d, -f4 | cut -d_ -f1 | sort | uniq > test.1
grep -v -f test.1 $folder.crams.run > test.2 
dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/BQSR/$folder/ | cut -d_ -f1 | sort | uniq > test.3
grep -f test.2 test.3 > bqsr.txt
for sample in $(cat bqsr.txt); do
   bam=$(dx find data --class file --path $KELLY:/CH_Exome/Workflow_Outputs/FINAL/BQSR/$folder/ --name "${sample}_23153_0_0.bqsr.bam" --brief | cut -d: -f2)

   dx run $KELLY:/test/tools/qc_exome_ukbb -iomni_vcf=file-G3zPG58J6XG2z5Qb8yGfPfP2 -ipicard_metric_accumulation_level="LIBRARY" -ibait_intervals=file-G5bgvB8JQ288QVGK2JP6BxYG -iper_base_intervals=file-G5bgvB8JQ288QVGK2JP6BxYG -itarget_intervals=file-G5bgvB8JQ288QVGK2JP6BxYG -iper_target_intervals=file-G5bgvB8JQ288QVGK2JP6BxYG -ireference=file-G3jKk3QJ6XG9KY8680qg9j8b -ireference_index=file-G3jKq0jJ6XGFqZ6Z9q4gB1yg -ibam=$bam -y --destination $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/$folder/ --name "${sample}-QC"
done



wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.*.*.tsv
########################################################################################################
########################################################################################################
########################################################################################################
"mem1_ssd1_v2_x8": 0.0528,
"mem2_ssd1_v2_x8": 0.0528,
"mem3_ssd1_v2_x8": 0.0528,
## run all since not many
for batch in {0000..0008}; do
   for file in a b; do
      word_count=$(wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.$batch.$file.tsv | awk '{print $1}')
      if [ $word_count -gt 1 ]; then
         dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.$batch.$file.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8
      else
      echo no
      fi
   done
done

analysis-G7kZV3QJ6XGGXgf9PZz3PKp3

"mem1_ssd1_v2_x8": 0.0528,
"mem2_ssd1_v2_x8": 0.0528, # default BQSR
"mem3_ssd1_v2_x8": 0.0528,
# 0000
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0000.a.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8 

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0000.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'

# 0001
sleep 3600
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0001.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 5=mem1_ssd1_v2_x8
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0001.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'

# 0002
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0002.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' 

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0002.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem2_ssd1_v2_x16 --instance-type 5=mem1_ssd1_v2_x8 
# 0003
sleep 7200
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0003.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 2=mem2_ssd1_v2_x4
sleep 1800
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0003.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 2=mem2_ssd1_v2_x8
# 0004

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0004.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'  --instance-type 5=mem2_ssd1_v2_x8
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0004.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 5=mem2_ssd1_v2_x4


# 0005

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0005.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0005.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem2_ssd1_v2_x8 --instance-type 5=mem2_ssd1_v2_x8


# 0006

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0006.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 4=mem2_ssd1_v2_x16


dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0006.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 5=mem2_ssd2_x8 --instance-type 4=mem2_ssd1_v2_x16


# 0007

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0007.a.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 5=mem2_ssd1_v2_x8

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0007.b.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 5=mem1_ssd1_v2_x8
# 0008
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0008.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'


cp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/48/batch.sh /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/50/batch.sh
cp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/48/dxworkflow.json /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/50/dxworkflow.json