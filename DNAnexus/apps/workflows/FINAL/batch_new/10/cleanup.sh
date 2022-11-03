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




folder=10
dx mkdir $NEW_EXOME:/Brian/Workflow_Outputs/MapQ0/$folder/Dups/
mkdir cleanup
(printf "State,Last modified,Size,Name,ID,\n";
dx ls -l --delimiter , $NEW_EXOME:/Brian/Workflow_Outputs/MapQ0/$folder/"*.vcf.gz*" | grep "[MK]B,${folder}") > cleanup/MapQ0.txt;

~/Bolton/UKBB/clean.up.crew.MapQ0.R --folder $folder

for file in $(cat cleanup/MapQ0.delete.txt); do 
      dx mv $file "$NEW_EXOME:/Brian/Workflow_Outputs/MapQ0/$folder/Dups/"
done



(printf "State,Last modified,Size,Name,ID,\n";
dx ls -l --delimiter , project-G8YxJ2QJP6z69ZVVFVfBYyPB:/Brian/Outputs/jie/"*.bam*" | grep "[MK]B,${folder}") > cleanup/complex.txt;

~/Bolton/UKBB/clean.up.crew.complex.R --folder none

for file in $(cat cleanup/complex.delete.txt); do 
      dx mv $file "project-G8YxJ2QJP6z69ZVVFVfBYyPB:/Brian/Outputs/jie/Dups/"
done