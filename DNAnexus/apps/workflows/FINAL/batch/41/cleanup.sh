folder=41
# cd /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder
# dx select $BRIAN

mkdir -p cleanup
# dx mkdir $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/$folder/Dups/
# dx mkdir $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/$folder/Dups/

(printf "State,Last modified,Size,Name,ID,\n";
dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/$folder/"*.txt";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Germline/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Lung/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Transplant/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Germline/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Lung/"*.txt" | grep "[MK]B,${folder}") > cleanup/Alignment_QC.txt;

for dir in MSK_Pileup MutectVEP VardictIntersectMutect; do 
   # dx mkdir $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/Dups/
   # dx mkdir $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/Dups/
   (printf "State,Last modified,Size,Name,ID,\n";
   dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/"*.vcf*";
   dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Germline/"*.vcf*" | grep "[MK]B,${folder}";
   dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Lung/"*.vcf*" | grep "[MK]B,${folder}";
   dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Transplant/"*.vcf*" | grep "[MK]B,${folder}";
   dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/Germline/"*.vcf*" | grep "[MK]B,${folder}";
   dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/Lung/"*.vcf*" | grep "[MK]B,${folder}") > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/cleanup/$dir.vcf.txt; 
done
~/Bolton/UKBB/clean.up.crew.R --folder $folder

wc -l cleanup/*.delete.txt
# mv brian
dx select $BRIAN
for dir in MSK_Pileup MutectVEP VardictIntersectMutect; do
   for file in $(cat cleanup/$dir.delete.txt); do 
      dx mv $file $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/Dups/
   done
   dx ls --brief $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/Dups/ > cleanup/brian.$dir.dupped
   grep -v -f cleanup/brian.$dir.dupped cleanup/$dir.delete.txt > cleanup/kelly.$dir.dupped
done
# dx select $BRIAN
# for dir in Alignment_QC; do
#    for file in $(cat cleanup/$dir.delete.txt); do 
#       dx mv $file $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/Dups/
#    done
#    dx ls --brief $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/Dups/ > cleanup/brian.$dir.dupped
#    grep -v -f cleanup/brian.$dir.dupped cleanup/$dir.delete.txt > cleanup/kelly.$dir.dupped
# done

dx select $KELLY
# chg proj and mv from kelly 
for dir in Alignment_QC MSK_Pileup MutectVEP VardictIntersectMutect; do
   for file in $(cat cleanup/kelly.$dir.dupped); do 
      dx mv $file $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/Dups/
   done
done
# for dir in Alignment_QC; do
#    for file in $(cat cleanup/kelly.$dir.dupped); do 
#       dx mv $file $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/Dups/
#    done
# done




############## Recheck there are 4*4010 QC, 1*4010 MSK, & 2*4010 Mutect/Vardict ##############
## 16040, 4010, 8020
wc -l $folder.crams.run


(printf "State,Last modified,Size,Name,ID,\n";
dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/$folder/"*.txt";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Germline/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Lung/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Transplant/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Germline/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Lung/"*.txt" | grep "[MK]B,${folder}") > cleanup/Alignment_QC.complete.txt;

for dir in MSK_Pileup MutectVEP VardictIntersectMutect; do 
   (printf "State,Last modified,Size,Name,ID,\n";
   dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/"*.vcf*";
   dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Germline/"*.vcf*" | grep "[MK]B,${folder}";
   dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Lung/"*.vcf*" | grep "[MK]B,${folder}";
   dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Transplant/"*.vcf*" | grep "[MK]B,${folder}";
   dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/Germline/"*.vcf*" | grep "[MK]B,${folder}";
   dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/Lung/"*.vcf*" | grep "[MK]B,${folder}") > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/cleanup/$dir.vcf.complete.txt; 
done
wc -l cleanup/*.complete*

