folder=13
cd /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder
for folder in 13; do
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

# for fold in {57..11}; do
#    dx mkdir project-G4qpk1jJQ285yvbXPFZKXkk8:/CH_Exome/Workflow_Outputs/FINAL/MSK_Pileup/$fold/Dups/
#    dx mv project-G4qpk1jJQ285yvbXPFZKXkk8:/CH_Exome/Workflow_Outputs/FINAL/MSK_Pileup/10/Dups/"$fold*.vcf.gz" project-G4qpk1jJQ285yvbXPFZKXkk8:/CH_Exome/Workflow_Outputs/FINAL/MSK_Pileup/$fold/Dups/
# done

# mv from kelly Alignment_QC
for folder in 10; do
   for dir in  MSK_Pileup MutectVEP VardictIntersectMutect; do
      for file in $(cat cleanup/$dir.delete.txt); do 
         dx mv $file "$NEW_EXOME:/Brian/Workflow_Outputs/$dir/$folder/Dups/"
      done
      dx ls --brief $NEW_EXOME:/Brian/Workflow_Outputs/$dir/$folder/Dups/ > cleanup/$dir.dupped
      grep -v -f cleanup/$dir.dupped cleanup/$dir.delete.txt
   done
done

wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/samples.new.txt
dx ls $NEW_EXOME:/Brian/Workflow_Outputs/$dir/$folder | cut -d_ -f1 | sort | uniq | wc -l
grep -f samples.new.txt samples.orginal.txt
dx ls 

############## Recheck there are 4x4002 QC, 1x4002 MSK, & 2x4002 Mutect/Vardict ##############
## 16009, 4003, 8005
wc -l $folder.crams.run

for folder in 10; do
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
      dx ls -l --delimiter , --full $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/"*.vcf*";
      dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Germline/"*.vcf*" | grep "[MK]B,${folder}";
      dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Lung/"*.vcf*" | grep "[MK]B,${folder}";
      dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Transplant/"*.vcf*" | grep "[MK]B,${folder}";
      dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/Germline/"*.vcf*" | grep "[MK]B,${folder}";
      dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/Lung/"*.vcf*" | grep "[MK]B,${folder}") > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder/cleanup/$dir.vcf.complete.txt; 
   done
done

wc -l cleanup/*.complete*
############## Recheck there are 4x4002 QC, 1x4002 MSK, & 2x4002 Mutect/Vardict ##############


########### compute1 ###########

folder=10
cd $EXOME/results/$folder/vcfs
for dir in MutectVEP VardictIntersectMutect; do 
   dx download $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/"*.vcf.gz";
   dx download $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Germline/"${folder}*.vcf.gz"
   dx download $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Lung/"${folder}*.vcf.gz"
   dx download $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Transplant/"${folder}*.vcf.gz"
   dx download $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/Germline/"${folder}*.vcf.gz"
   dx download $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/Lung/"${folder}*.vcf.gz"
done




for fold in {10..60}; do
   cd $EXOME/results/$fold
   rm *.tsv; rm *.log
   cp ~/nSamples.sh .
   chmod u+x nSamples.sh
   cp ~/combine.sh .
   chmod u+x combine.sh
   mkdir -p vcfs
   mkdir -p final_outputs
   mkdir -p combined
done

