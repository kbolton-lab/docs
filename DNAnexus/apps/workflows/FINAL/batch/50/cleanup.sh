folder=50
mkdir -p cleanup

(printf "State,Last modified,Size,Name,ID,\n";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/$folder/"*.txt";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Germline/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Lung/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Transplant/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Germline/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Lung/"*.txt" | grep "[MK]B,${folder}") > cleanup/Alignment_QC.txt;

for dir in MSK_Pileup MutectVEP VardictIntersectMutect; do 
   (printf "State,Last modified,Size,Name,ID,\n";
   dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/"*.vcf*";
   dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Germline/"*.vcf*" | grep "[MK]B,${folder}";
   dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Lung/"*.vcf*" | grep "[MK]B,${folder}";
   dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Transplant/"*.vcf*" | grep "[MK]B,${folder}";
   dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/Germline/"*.vcf*" | grep "[MK]B,${folder}";
   dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/Lung/"*.vcf*" | grep "[MK]B,${folder}") > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/cleanup/$dir.vcf.txt; 
done
~/Bolton/UKBB/clean.up.crew.R --folder $folder
wc -l cleanup/*.delete.txt

# mv kelly
dx select $KELLY
for dir in Alignment_QC MSK_Pileup MutectVEP VardictIntersectMutect; do
   for file in $(cat cleanup/$dir.delete.txt); do 
      dx mv $file "$KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/Dups/"
   done
   dx ls --brief $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/Dups/ > cleanup/kelly.$dir.dupped
   grep -v -f cleanup/kelly.$dir.dupped cleanup/$dir.delete.txt > cleanup/brian.$dir.dupped
done
# for dir in Alignment_QC; do
#    for file in $(cat cleanup/$dir.delete.txt); do 
#       dx mv $file "$KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/Dups/"
#    done
#    dx ls --brief $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/Dups/ > cleanup/kelly.$dir.dupped
#    grep -v -f cleanup/kelly.$dir.dupped cleanup/$dir.delete.txt > cleanup/brian.$dir.dupped
# done


dx select $BRIAN
# chg proj and mv from brian 
for dir in MSK_Pileup MutectVEP VardictIntersectMutect; do
   for file in $(cat cleanup/brian.$dir.dupped); do 
      dx mv $file "$BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/Dups/"
   done
done
# for dir in Alignment_QC; do
#    for file in $(cat cleanup/brian.$dir.dupped); do 
#       dx mv $file "$BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/Dups/"
#    done
# done



############## Recheck there are 4*3909 QC, 1*3909 MSK, & 2*3909 Mutect/Vardict ##############
## 15736, 3909, 7818
wc -l $folder.crams.run
cd /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder

(printf "State,Last modified,Size,Name,ID,\n";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/$folder/"*.txt";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Germline/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Lung/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Transplant/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Germline/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Lung/"*.txt" | grep "[MK]B,${folder}") > cleanup/Alignment_QC.complete.txt;

for dir in MSK_Pileup MutectVEP VardictIntersectMutect; do 
   (printf "State,Last modified,Size,Name,ID,\n";
   dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/"*.vcf*";
   dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Germline/"*.vcf*" | grep "[MK]B,${folder}";
   dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Lung/"*.vcf*" | grep "[MK]B,${folder}";
   dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Transplant/"*.vcf*" | grep "[MK]B,${folder}";
   dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/Germline/"*.vcf*" | grep "[MK]B,${folder}";
   dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/Lung/"*.vcf*" | grep "[MK]B,${folder}") > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/cleanup/$dir.vcf.complete.txt; 
done
wc -l cleanup/*.complete*

