function TOTALK {
   PROJECT=$KELLY
   folder=$1

   (dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/MSK_Pileup/"Lung"/"*.normal.pileup.*" | grep ^$folder; 
    dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/MSK_Pileup/"Germline"/"*.normal.pileup.*" | grep ^$folder;
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MSK_Pileup/"Lung"/"*.normal.pileup.*" | grep ^$folder; 
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MSK_Pileup/"Germline"/"*.normal.pileup.*" | grep ^$folder; 
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MSK_Pileup/"Transplant"/"*.normal.pileup.*" | grep ^$folder;
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MSK_Pileup/$folder/"*.normal.pileup.*") | cut -d_ -f1 | sort | uniq -c | sort -k1,1rn | awk '$1 > 1' > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dups.txt

   (dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"Lung"/"*.vcf.gz"; 
    dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"Germline"/"*.vcf.gz") | cut -d_ -f1 | sort | uniq | grep ^$folder > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/brian.vardict.LG.total.txt

   (dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"Lung"/"*.vcf.gz"; 
    dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"Germline"/"*.vcf.gz") | cut -d_ -f1 | sort | uniq | grep ^$folder > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/brian.mutect.LG.total.txt
   
   (dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"Lung"/"*.vcf.gz"; 
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"Germline"/"*.vcf.gz"; 
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"Transplant"/"*.vcf.gz") | cut -d_ -f1 | sort | uniq | grep ^$folder > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/kelly.vardict.TLG.total.txt


   dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"$folder"/"*.vcf.gz" | cut -d_ -f1 | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/kelly.vardict.total.txt

   (dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"Lung"/"*.vcf.gz"; 
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"Germline"/"*.vcf.gz"; 
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"Transplant"/"*.vcf.gz") | cut -d_ -f1 | sort | uniq | grep ^$folder > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/kelly.mutect.TLG.total.txt
 
   dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"$folder"/"*.vcf.gz" | cut -d_ -f1 | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/kelly.mutect.total.txt

   cat /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/*.vardict*.total.txt | grep ^$folder | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/VARDICT.total.txt

   cat /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/*.mutect*.total.txt | grep ^$folder | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/MUTECT.total.txt

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
   echo "\n"
   wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.*.*.tsv
   ## DONT DUPLICATE!!!!
   cat /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.*.*.tsv | grep -v batch | awk '{print $1}' > resend.txt
   grep -f resend.txt running.txt
   
   total_to_run=$(wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.*.*.tsv | tail -n1 | awk -v remove=$i '{print $1-remove}')
   total_ran_run=$(wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/total.RAN.RUNNING.txt | awk '{print $1}')
   total_crams=$(wc -l  /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/$folder.crams.run | awk '{print $1}')
   echo "$total_to_run to run + $total_ran_run ran = $(echo $total_to_run + $total_ran_run | bc) total"
   echo "total $total_crams crams"
}

function TOTALB {
   folder=$1
   PROJECT=$BRIAN

   (dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/MSK_Pileup/"Lung"/"*.normal.pileup.*" | grep ^$folder; 
    dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/MSK_Pileup/"Germline"/"*.normal.pileup.*" | grep ^$folder;
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MSK_Pileup/"Lung"/"*.normal.pileup.*" | grep ^$folder; 
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MSK_Pileup/"Germline"/"*.normal.pileup.*" | grep ^$folder; 
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MSK_Pileup/"Transplant"/"*.normal.pileup.*" | grep ^$folder;
    dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/MSK_Pileup/$folder/"*.normal.pileup.*") | cut -d_ -f1 | sort | uniq -c | sort -k1,1rn | awk '$1 > 1' > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dups.txt

   (dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"Lung"/"*.vcf.gz"; 
    dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"Germline"/"*.vcf.gz") | cut -d_ -f1 | sort | uniq | grep ^$folder > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/brian.vardict.LG.total.txt

   dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"$folder"/"*.vcf.gz" | cut -d_ -f1 | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/brian.vardict.total.txt

   (dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"Lung"/"*.vcf.gz"; 
    dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"Germline"/"*.vcf.gz") | cut -d_ -f1 | sort | uniq | grep ^$folder > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/brian.mutect.LG.total.txt

   dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"$folder"/"*.vcf.gz" | cut -d_ -f1 | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/brian.mutect.total.txt


   (dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"Lung"/"*.vcf.gz"; 
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"Germline"/"*.vcf.gz"; 
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"Transplant"/"*.vcf.gz") | cut -d_ -f1 | sort | uniq | grep ^$folder > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/kelly.vardict.TLG.total.txt

   (dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"Lung"/"*.vcf.gz"; 
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"Germline"/"*.vcf.gz"; 
    dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"Transplant"/"*.vcf.gz") | cut -d_ -f1 | sort | uniq | grep ^$folder > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/kelly.mutect.TLG.total.txt
 
   
   cat /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/*.vardict*.total.txt | grep ^$folder | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/VARDICT.total.txt

   cat /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/*.mutect*.total.txt | grep ^$folder | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/MUTECT.total.txt

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
   echo "\n"
   wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.*.*.tsv
   ## DONT DUPLICATE!!!!
   cat /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.*.*.tsv | grep -v batch | awk '{print $1}' > resend.txt
   grep -f resend.txt running.txt

   total_to_run=$(wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/dx_batch.*.*.tsv | tail -n1 | awk -v remove=$i '{print $1-remove}')
   total_ran_run=$(wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/total.RAN.RUNNING.txt | awk '{print $1}')
   total_crams=$(wc -l  /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/$folder.crams.run | awk '{print $1}')
   echo "$total_to_run to run + $total_ran_run ran = $(echo $total_to_run + $total_ran_run | bc) total"
   echo "total $total_crams crams"
}