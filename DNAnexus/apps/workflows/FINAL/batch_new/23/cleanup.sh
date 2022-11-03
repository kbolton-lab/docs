folder=23
cd /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder
for folder in 23; do
   mkdir -p cleanup
   dx mkdir $NEW_EXOME:/Brian/Workflow_Outputs/Alignment_QC/$folder/Dups/
   (printf "State,Last modified,Size,Name,ID,\n";
   dx ls -l --delimiter , $NEW_EXOME:/Brian/Workflow_Outputs/Alignment_QC/$folder/"*.txt" | grep "[MK]B,${folder}") > cleanup/Alignment_QC.txt;

   for dir in MSK_Pileup MutectVEP VardictIntersectMutect; do 
      
      dx mkdir $NEW_EXOME:/Brian/Workflow_Outputs/$dir/$folder/Dups/
      (printf "State,Last modified,Size,Name,ID,\n";
      dx ls -l --delimiter , $NEW_EXOME:/Brian/Workflow_Outputs/$dir/$folder/"*.vcf*" | grep "[MK]B,${folder}") > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/cleanup/$dir.vcf.txt; 
   done
   ~/Bolton/UKBB/clean.up.crew.R --folder $folder
done


for dir in Alignment_QC MSK_Pileup MutectVEP VardictIntersectMutect; do
   for file in $(cat cleanup/$dir.delete.txt); do 
      dx mv $file "$NEW_EXOME:/Brian/Workflow_Outputs/$dir/$folder/Dups/"
   done
   dx ls --brief $NEW_EXOME:/Brian/Workflow_Outputs/$dir/$folder/Dups/ > cleanup/$dir.dupped
   grep -v -f cleanup/$dir.dupped cleanup/$dir.delete.txt
done


wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/samples.new.txt
dx ls $NEW_EXOME:/Brian/Workflow_Outputs/$dir/$folder | cut -d_ -f1 | sort | uniq | wc -l
grep -f samples.new.txt samples.orginal.txt
dx ls 



