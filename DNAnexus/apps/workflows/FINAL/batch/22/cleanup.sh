folder=22
# cd /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/$folder
# dx select $KELLY

mkdir -p cleanup
# dx mkdir $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/$folder/Dups/
# dx mkdir $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/$folder/Dups/

(printf "State,Last modified,Size,Name,ID,\n";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/$folder/"*.txt";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Germline/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Lung/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Transplant/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Germline/"*.txt" | grep "[MK]B,${folder}";
dx ls -l --delimiter , $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/Lung/"*.txt" | grep "[MK]B,${folder}") > cleanup/Alignment_QC.txt;

for dir in MSK_Pileup MutectVEP VardictIntersectMutect; do 
   # dx mkdir $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/Dups/
   # dx mkdir $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/Dups/
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
tail -n+2 cleanup/Alignment_QC.txt | cut -d, -f4 | cut -d_ -f1 | sort | uniq > test.2
dx ls project-G4qpk1jJQ285yvbXPFZKXkk8:/CH_Exome/Workflow_Outputs/FINAL/BQSR/18/ | cut -d_ -f1 | sort | uniq > test.3
# mv kelly
dx select $KELLY
for dir in  MSK_Pileup MutectVEP VardictIntersectMutect; do
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

wc -l cleanup/brian.*.dupped
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
## 15736, 3966, 7932
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

############## Recheck there are 4x4002 QC, 1x4002 MSK, & 2x4002 Mutect/Vardict ##############


########### compute1 ###########

folder=16
cd $EXOME/results/$folder/vcfs
for dir in MutectVEP VardictIntersectMutect; do 
   dx download $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/"*.vcf.gz";
   dx download $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Germline/"${folder}*.vcf.gz"
   dx download $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Lung/"${folder}*.vcf.gz"
   dx download $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Transplant/"${folder}*.vcf.gz"
   dx download $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/Germline/"${folder}*.vcf.gz"
   dx download $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/Lung/"${folder}*.vcf.gz"
done
bsub -oo mutect.log -G compute-bolton -q general -M 4G -R 'rusage[mem=4G]' -a 'docker(kboltonlab/dnanexus:1.0)' /bin/bash ./download.sh 16 MutectVEP
bsub -oo vardict.log -G compute-bolton -q general -M 4G -R 'rusage[mem=4G]' -a 'docker(kboltonlab/dnanexus:1.0)' /bin/bash ./download.sh 16 VardictIntersectMutect
#!/bin/bash
# download.sh
DX2
folder=$1
dir=$2
dx download $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/"*.vcf.gz" -f
dx download $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Germline/"${folder}*.vcf.gz" -f
dx download $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Lung/"${folder}*.vcf.gz" -f
dx download $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$dir/Transplant/"${folder}*.vcf.gz" -f
dx download $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/Germline/"${folder}*.vcf.gz" -f
dx download $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/Lung/"${folder}*.vcf.gz" -f




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

