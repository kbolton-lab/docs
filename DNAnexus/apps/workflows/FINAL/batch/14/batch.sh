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

folder=14
PROJECT=$KELLY
workflow=workflow-G7bjpv0JQ2899YVx8q6PZYFZ
## 2 stage to get vardict
workflow=workflow-G7jybvQJQ288G8FKJGz7xG1q
dx generate_batch_inputs -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam='(.*)_23153_0_0.cram$' -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index='(.*)_23153_0_0.cram.crai' --path "UKBB_Exome_2021:Bulk/Exome sequences/Exome OQFE CRAM files/$folder" -o "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch"

# "mem1_ssd1_v2_x8": 0.0528,
# "mem2_ssd1_v2_x8": 0.0528,
# "mem3_ssd1_v2_x8": 0.0528,



# # BQSR 
# # tp_53_bqsr
# dx ls project-G4qpk1jJQ285yvbXPFZKXkk8:/CH_Exome/Workflow_Outputs/FINAL/BQSR/Germline_TP53/"*.bam" | cut -d_ -f1 | grep ^$folder | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/tp_53_bqsr.txt

# dx ls project-G4qpk1jJQ285yvbXPFZKXkk8:/CH_Exome/Workflow_Outputs/FINAL/BQSR/$folder/"*.bam" | cut -d_ -f1 | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/bqsr.txt

# cat /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.*.*.tsv | grep -v batch | awk '{print $1}' > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/ids_to_run.txt

# # for bam in $(grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/ids_to_run.txt /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/tp_53_bqsr.txt); do 
# #    dx rm project-G4qpk1jJQ285yvbXPFZKXkk8:/CH_Exome/Workflow_Outputs/FINAL/BQSR/Germline_TP53/"$bam*.bam*"
# # done

# for bam in $(grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/ids_to_run.txt /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/bqsr.txt); do 
#    dx rm project-G4qpk1jJQ285yvbXPFZKXkk8:/CH_Exome/Workflow_Outputs/FINAL/BQSR/$folder/"$bam*.bam*"
# done


folder=14
workflow=workflow-G7jybvQJQ288G8FKJGz7xG1q
rm VIM.app.run
touch VIM.app.run
wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/14/vardict.intersect.mutect.need.to.run.txt
for sample in $(cat /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/14/vardict.intersect.mutect.need.to.run.txt2); do
   mutect_file=$(dx find data --class file --path $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/ --name "${sample}_23153_0_0.mutect.vep.annotated.vcf.gz" --brief | cut -d: -f2)
   mutect_index=$(dx find data --class file --path $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/ --name "${sample}_23153_0_0.mutect.vep.annotated.vcf.gz.tbi"  --brief | cut -d: -f2)
   vardict_file=$(dx find data --class file --path $KELLY:/CH_Exome/Workflow_Outputs/FINAL/VardictJava/ --name "${sample}_23153_0_0.vardict.BCBIOfiltered.vcf.gz" --brief | cut -d: -f2)
   vardict_index=$(dx find data --class file --path $KELLY:/CH_Exome/Workflow_Outputs/FINAL/VardictJava/ --name "${sample}_23153_0_0.vardict.BCBIOfiltered.vcf.gz.tbi" --brief | cut -d: -f2)

   dx run $KELLY:/test/tools/vardict_pon2at2percent_final2 -iintersect_vcf=$mutect_file -iintersect_vcf_index=$mutect_index -ivcf2PON=file-G7B1gyQJ6XG85b9j1k11k0FQ -ivcf2PON_index=file-G7B1j68J6XG6z5ZF24vvk2qX -ivcf=$vardict_file -ivcf_index=$vardict_index -ireference=file-G3jKk3QJ6XG9KY8680qg9j8b -ireference_index=file-G3jKq0jJ6XGFqZ6Z9q4gB1yg -icaller=vardict -y --destination $KELLY:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/14/ --name "${sample}-vardict_pon2at2percent_final2"
done

########################################################################################################
########################################################################################################
########################################################################################################
"mem1_ssd1_v2_x8": 0.0528,
"mem2_ssd1_v2_x8": 0.0528,
"mem3_ssd1_v2_x8": 0.0528,
## run all since not many
for batch in {0000..0007}; do
   for file in a b; do
      word_count=$(wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.$batch.$file.tsv | awk '{print $1}')
      if [ $word_count -gt 1 ]; then
         dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.$batch.$file.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem2_ssd1_v2_x8 --instance-type 5=mem2_ssd1_v2_x8
      else
      echo no
      fi
   done
done

"mem1_ssd1_v2_x8": 0.0528,
"mem2_ssd1_v2_x8": 0.0528, # default BQSR
"mem3_ssd1_v2_x8": 0.0528,
# 0000
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0000.a.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem1_ssd1_v2_x16 --instance-type 8=mem1_ssd1_v2_x8
sleep 1200; dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0000.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'

# 0001
sleep 600; dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0001.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 8=mem1_ssd1_v2_x8 --instance-type 5=mem1_ssd1_v2_x8
sleep 2400; dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0001.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem1_ssd1_v2_x16 

# 0002
sleep 2400; 
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0002.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 8=mem1_ssd1_v2_x8
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0002.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem2_ssd1_v2_x16 --instance-type 5=mem1_ssd1_v2_x8 --instance-type 8=mem1_ssd1_v2_x8 
# 0003
sleep 2400; 
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0003.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0003.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_hdd2_v2_x4 --instance-type 8=mem1_ssd1_v2_x8
# 0004
sleep 2400
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0004.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 8=mem2_ssd2_x8 --instance-type 5=mem2_ssd1_v2_x8
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0004.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_hdd2_v2_x4 --instance-type 5=mem2_ssd1_v2_x4


# 0005

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0005.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 8=mem1_ssd1_v2_x8
sleep 2400
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0005.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem2_ssd1_v2_x8 --instance-type 5=mem2_ssd1_v2_x8


# 0006
sleep 2100
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0006.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0006.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 5=mem2_ssd2_x8


# 0007
sleep 600
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0007.a.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 5=mem2_ssd1_v2_x8

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0007.b.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 8=mem1_ssd1_v2_x8
# 0008
# sleep 600; dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0008.tsv   -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'


cp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/20/batch.sh /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/22/batch.sh
cp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/20/dxworkflow.json /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/22/dxworkflow.json



