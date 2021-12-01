dx generate_batch_inputs -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam='(.*)_23153_0_0.cram$' -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index='(.*)_23153_0_0.cram.crai' --path "UKBB_Exome_2021:Bulk/Exome sequences/Exome OQFE CRAM files/12" -o "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/dx_batch"

(dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Lung/"*.final.tsv" | cut -d_ -f1 | sort | uniq | grep "^12";
dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Germline/"*.final.tsv" | cut -d_ -f1 | sort | uniq | grep "^12";
dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/12/"*.final.tsv" | cut -d_ -f1 | sort | uniq | grep "^12") > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/brian.12.txt

(dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Lung/"*.final.tsv" | cut -d_ -f1 | sort | uniq | grep "^12";
dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Germline/"*.final.tsv" | cut -d_ -f1 | sort | uniq | grep "^12";
dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Transplant/"*.final.tsv" | cut -d_ -f1 | sort | uniq | grep "^12";
dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/12/"*.final.tsv" | cut -d_ -f1 | sort | uniq | grep "^12") > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/kelly.12.txt


wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/dx_batch.000*.tsv

for file in $(ls /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/dx_batch.000*.tsv | grep -v remaining | grep -v ran); do
    name=$(basename $file .tsv)
    echo $name
    grep -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/kelly.12.txt $file > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/$name.ran.tsv
done

wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/dx_batch.000*.remaining.tsv

KELLY
/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/dx_batch.0003.remaining.tsv
/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/dx_batch.0005.remaining.tsv
/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/dx_batch.0007.remaining.tsv


dx run workflow-G6Q7gZ0JQ281kBv3GQxFGf5g --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/dx_batch.0003.remaining.tsv -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$KELLY -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
dx run workflow-G6Q7gZ0JQ281kBv3GQxFGf5g --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/dx_batch.0005.remaining.tsv -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$KELLY -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
dx run workflow-G6Q7gZ0JQ281kBv3GQxFGf5g --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/dx_batch.0007.remaining.tsv -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$KELLY -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'


BRIAN
/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/dx_batch.0004.remaining.tsv
/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/dx_batch.0006.remaining.tsv

workflow-G6Q7gZ0JQ281kBv3GQxFGf5g

dx run workflow-G6Q7gZ0JQ281kBv3GQxFGf5g --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/dx_batch.0004.remaining.tsv -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$BRIAN -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'
dx run workflow-G6Q7gZ0JQ281kBv3GQxFGf5g --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/12/FINAL_COMBINE_JOBS_12/dx_batch.0006.remaining.tsv -istage-G4KqPz0J6XGJZGB842qJVYQK.project=$BRIAN -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}'



(dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Lung/"*.final.tsv" | cut -d_ -f1 | sort | uniq;
dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Germline/"*.final.tsv" | cut -d_ -f1 | sort | uniq;
dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/12/"*.final.tsv" | cut -d_ -f1 | sort | uniq;
dx ls $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/13/"*.final.tsv" | cut -d_ -f1 | sort | uniq) > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/brian.total.txt

(dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Lung/"*.final.tsv" | cut -d_ -f1 | sort | uniq;
dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Germline/"*.final.tsv" | cut -d_ -f1 | sort | uniq;
dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/Transplant/"*.final.tsv" | cut -d_ -f1 | sort | uniq;
dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/12/"*.final.tsv" | cut -d_ -f1 | sort | uniq | grep "^12") > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/kelly.total.txt

wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/*.total.txt