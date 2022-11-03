

      (printf "State,Last modified,Size,Name,ID,\n";
      dx ls -l --delimiter , project-G3Yj1vjJ6XG579jbKyjXPGGY:/CH_Exome/Workflow_Outputs/FINAL/Mutect_Rerun/Brian/Workflow_Outputs/FilterMutectCalls/$folder/"*.vcf*" | grep "[MK]B,${folder}") > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dups.vcf.txt; 

   ~/Bolton/UKBB/clean.up.crew.R --folder $folder
   ~/Bolton/UKBB/clean.up.crew3.R --folder $folder

for file in $(cat cleanup.delete.txt); do 
      dx mv $file "project-G3Yj1vjJ6XG579jbKyjXPGGY:/CH_Exome/Workflow_Outputs/FINAL/Mutect_Rerun/Brian/Workflow_Outputs/FilterMutectCalls/$folder/Dups/"
done






