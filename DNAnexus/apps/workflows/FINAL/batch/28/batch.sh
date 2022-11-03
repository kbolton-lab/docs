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

folder=28
PROJECT=$KELLY
workflow=workflow-G7b8BPQJQ288V8Q435J96151
dx generate_batch_inputs -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam='(.*)_23153_0_0.cram$' -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index='(.*)_23153_0_0.cram.crai' --path "UKBB_Exome_2021:Bulk/Exome sequences/Exome OQFE CRAM files/$folder" -o "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch"

# "mem1_ssd1_v2_x8": 0.0528,
# "mem2_ssd1_v2_x8": 0.0528,
# "mem3_ssd1_v2_x8": 0.0528,

function TOTAL {
   folder=$1

   (dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"Lung"/"*.vcf.gz"; 
    dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"Germline"/"*.vcf.gz") | cut -d_ -f1 | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/brian.vardict.total.txt

   (dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"Lung"/"*.vcf.gz"; 
    dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"Germline"/"*.vcf.gz") | cut -d_ -f1 | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/brian.mutect.total.txt


   (dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"$folder"/"*.vcf.gz";
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"Lung"/"*.vcf.gz"; 
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"Germline"/"*.vcf.gz"; 
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"Transplant"/"*.vcf.gz") | cut -d_ -f1 | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/kelly.vardict.total.txt
 
   (dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"$folder"/"*.vcf.gz";
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"Lung"/"*.vcf.gz"; 
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"Germline"/"*.vcf.gz"; 
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"Transplant"/"*.vcf.gz") | cut -d_ -f1 | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/kelly.mutect.total.txt

   cat /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/*.vardict.total.txt | grep ^$folder | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/VARDICT.total.txt

   cat /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/*.mutect.total.txt | grep ^$folder | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/MUTECT.total.txt

   # intersect to see which are done # 853 either order
   grep -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/VARDICT.total.txt /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/MUTECT.total.txt > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/total.RAN.txt

   # not done/failed?
   # grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/MUTECT.total.txt /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/VARDICT.total.txt

   # Any in progress?
   # https://askubuntu.com/questions/642563/passing-variable-into-awk-gsub
   dx find executions --project $PROJECT --state in_progress -n 2000 --origin-jobs --name "$folder-*" | grep analysis | awk -v folder=$folder '{gsub("'"$folder-"'","",$2); print $2}' > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/running.txt
   
   cat /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/total.RAN.txt /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/running.txt | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/total.RAN.RUNNING.txt
   # /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/workflow.total.txt

   i=0
   for batch in {0000..0008}; do
      head -n251 /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.$batch.tsv | grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/total.RAN.RUNNING.txt > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.$batch.a.tsv
      wc_l=$(wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.$batch.tsv | awk '{print $1}')
      if [[ $wc_l -gt 1 ]]; then
         echo "increment"
         i=$((i+1))
      else 
         echo "no increment"
      fi
      
      if [[ $wc_l -gt 251 ]]; then
         echo "do this"
         i=$((i+1))
         my_tail=$(echo $wc_l-250-1 | bc)
         (head -n1 /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.$batch.tsv; tail -n${my_tail} /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.$batch.tsv | grep -v -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/total.RAN.RUNNING.txt) > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.$batch.b.tsv
      else 
         echo "don't do this"
      fi
      wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.$batch.*.tsv
   done
   wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.*.*.tsv

   total_to_run=$(wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.*.*.tsv | tail -n1 | awk -v remove=$i '{print $1-remove}')
   total_ran_run=$(wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/total.RAN.RUNNING.txt | awk '{print $1}')
   total_crams=$(wc -l  /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/$folder.crams.run | awk '{print $1}')
   echo "$total_to_run to run + $total_ran_run ran = $(echo $total_to_run + $total_ran_run | bc) total"
   echo "total $total_crams crams"
}
export -f TOTAL

TOTAL 28


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
         dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.$batch.$file.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 5=mem2_ssd2_x8
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
sleep 1200; dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0000.b.tsv  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$PROJECT -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'

# 0001
sleep 600; dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0001.a.tsv  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$PROJECT -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 8=mem1_ssd1_v2_x8 --instance-type 5=mem1_ssd1_v2_x8
sleep 2400; dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0001.b.tsv  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$PROJECT -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem1_ssd1_v2_x16 

# 0002
sleep 2400; 
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0002.a.tsv  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$PROJECT -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 8=mem1_ssd1_v2_x8
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0002.b.tsv  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$PROJECT -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem2_ssd1_v2_x16 --instance-type 5=mem1_ssd1_v2_x8 --instance-type 8=mem1_ssd1_v2_x8 
# 0003
sleep 2400; 
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0003.a.tsv  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$PROJECT -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0003.b.tsv  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$PROJECT -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_hdd2_v2_x4 --instance-type 8=mem1_ssd1_v2_x8
# 0004
sleep 2400
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0004.a.tsv  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$PROJECT -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 8=mem2_ssd2_x8 --instance-type 5=mem2_ssd1_v2_x8
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0004.b.tsv  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$PROJECT -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_hdd2_v2_x4 --instance-type 5=mem2_ssd1_v2_x4


# 0005

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0005.a.tsv  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$PROJECT -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8  --instance-type 8=mem1_ssd1_v2_x8
sleep 2400
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0005.b.tsv  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$PROJECT -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem2_ssd1_v2_x8 --instance-type 5=mem2_ssd1_v2_x8


# 0006
sleep 2100
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0006.a.tsv  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$PROJECT -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0006.b.tsv  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$PROJECT -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 5=mem2_ssd2_x8


# 0007
sleep 600
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0007.a.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 5=mem2_ssd1_v2_x8

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0007.b.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --instance-type 0=mem3_ssd1_v2_x8 --instance-type 5=mem1_ssd1_v2_x8
# 0008
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.0008.tsv  -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$PROJECT -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'


cp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/20/batch.sh /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/22/batch.sh
cp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/20/dxworkflow.json /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/22/dxworkflow.json