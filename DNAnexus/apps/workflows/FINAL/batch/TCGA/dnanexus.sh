workflow=workflow-G7jgVZ0Jf313xQ12J94gZ7Jp # old wrong name, still use initially
# workflow= # new correct name, use later
# workflow=workflow-G7vkJ0QJf315kfqG6FBZ8x79 # increase bqsr timeout
# workflow=workflow-G7vp718Jf31PJb35569p16ZV # increase bqsr timeout 
PROJECT=project-G428g40Jf312vx5727kY9KXZ

for uuid in $(grep -v -f blood.bams.downloaded 500.txt | awk '{print $1}'); do dx run $JIE:/test/tools/tcga_dl -y -iuuid=$uuid --destination project-G428g40Jf312vx5727kY9KXZ:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/ --name "tcga_dl-${uuid}"; done


dx generate_batch_inputs -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam='(.*).bam$' -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index='(.*).bai$' --path "project-G428g40Jf312vx5727kY9KXZ:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test" -o "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/TCGA/batch/dx_batch_blood"
dx generate_batch_inputs -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam='(.*).bam$' -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index='(.*).bai$' --path "project-G428g40Jf312vx5727kY9KXZ:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/multiples/merged" -o "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/TCGA/batch/dx_batch_blood_merged"

dx generate_batch_inputs -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam='(.*).bam$' -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index='(.*).bai$' --path "project-G428g40Jf312vx5727kY9KXZ:/CH_Exome/Inputs/TCGA/BAMS/Tumor_Test" -o "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/TCGA/dx_batch_tumor"

BQSR_delete
dx ls -l project-G428g40Jf312vx5727kY9KXZ:/CH_Exome/Workflow_Outputs/FINAL/BQSR/ | awk '{gsub(/\(|\)/,"",$7); print $6,$7}' OFS='\t' | tail -n+4 > bqsr.txt
cat blood_done tumor_done > all.done 
grep -f all.done bqsr.txt | while read file id; do
   echo $file
   dx mv $id $JIE:/CH_Exome/Workflow_Outputs/FINAL/BQSR/Remove/
done

721-77=644
1337 D-DNA bams 
1644-1337
for sample in $(cat blood_done); do
   qc=$(grep ${sample} qc.txt | wc -l)
   printf "$sample\t$qc\n"
done


C509.TCGA-78-7160-10A-01D-2036-08.1_gdc_realn
C508.TCGA-77-8154-10A-01D-2244-08.1_gdc_realn
# Normal Blood
# /[A-Z0-9]{4}-10[A-Z]
dx ls project-G428g40Jf312vx5727kY9KXZ:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/"*.bam" | grep -E "[A-Z0-9]{4}-10[A-Z]-[0-9]{2}D" | cut -d- -f1,2,3,4 > blood_D-DNA_bams
for uuid in $(grep -v -f blood_D-DNA_bams blood_normal.D_DNA.tsv | awk '{print $1}'); do dx run $JIE:/test/tools/tcga_dl -y -iuuid=$uuid --destination $JIE:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/ --name tcga_dl-${uuid}; sleep 25; done

#!/usr/bin/bash
source /usr/local/dx-toolkit/environment;
#dx ls project-G428g40Jf312vx5727kY9KXZ:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/"*.bam" | grep -E "[A-Z0-9]{4}-10[A-Z]-[0-9]{2}D" | cut -d- -f1,2,3,4 > blood_D-DNA_bams
(dx ls project-G428g40Jf312vx5727kY9KXZ:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/"*.bam"; dx ls project-G428g40Jf312vx5727kY9KXZ:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/multiples/"*.bam") > blood_D-DNA_bams
for uuid in $(grep -v -f blood_D-DNA_bams blood_normal.D_DNA.multiple.tsv | awk '{print $1}'); do dx run project-G428g40Jf312vx5727kY9KXZ:/test/tools/tcga_dl -y -iuuid=$uuid --destination project-G428g40Jf312vx5727kY9KXZ:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/multiples/ --name tcga_dl-${uuid}; sleep 15; done

bsub -oo dna.log -G compute-bolton -g /jietest -q general -M 32G -R 'rusage[mem=32G]' -a 'docker(kboltonlab/dnanexus:1.0)' ./script.single.sh
bsub -oo dna1.log -G compute-bolton -g /jietest -q general -M 32G -R 'rusage[mem=32G]' -a 'docker(kboltonlab/dnanexus:1.0)' ./script.single.noname.sh
bsub -oo dna2.log -G compute-bolton -g /jietest -q general -M 32G -R 'rusage[mem=32G]' -a 'docker(kboltonlab/dnanexus:1.0)' ./script.multiple.sh
bsub -oo dna3.log -G compute-bolton -g /jietest -q general -M 32G -R 'rusage[mem=32G]' -a 'docker(kboltonlab/dnanexus:1.0)' ./script.multiple.noname.sh
/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/
/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/multiples/

tail -n+11 blood_normal.D_DNA.no.name.tsv | awk -F'\t' '{print $1,$2,$3}' | 
paste blood_normal.D_DNA.no.name.multiple.tsv dug.t | awk '{print $2,$3,$9}' | while read file_name sample number; do
   file_prefix=$(basename $file_name .bam)
   dx ls $JIE:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/multiples/$file_prefix.bam
   dx ls $JIE:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/multiples/$file_prefix.bai
   dx mv $JIE:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/multiples/$file_prefix.bam $JIE:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/multiples/${sample}_${number}_gdc_realn.bam
   dx mv $JIE:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/multiples/$file_prefix.bai $JIE:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/multiples/${sample}_gdc_realn.bai
   echo 
done
dx ls project-G428g40Jf312vx5727kY9KXZ:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/multiples/"*.bam" > multiples.txt

cat multiples.2.txt | while read bam_in1 bam_in2; do
   sample=$(echo $bam_in1 | cut -d_ -f1);
   prefix1=$(basename $bam_in1 .bam)
   prefix2=$(basename $bam_in2 .bam)
   bam1=$(dx find data --class file --path $JIE:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/multiples/ --name "$bam_in1" --brief | cut -d: -f2)
   bam2=$(dx find data --class file --path $JIE:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/multiples/ --name "$bam_in2" --brief | cut -d: -f2)
   bai1=$(dx find data --class file --path $JIE:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/multiples/ --name "$prefix1.bai" --brief | cut -d: -f2)
   bai2=$(dx find data --class file --path $JIE:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/multiples/ --name "$prefix2.bai" --brief | cut -d: -f2)

   dx run $JIE:/test/tools/merge_bam -ibam_in1=$bam1 -ibam_in2=$bam2 -ibai_in1=$bai1 -ibai_in2=$bai2 -iout_name=${sample}_merged_gdc_realn --destination $JIE:/CH_Exome/Inputs/TCGA/BAMS/Blood_Test/multiples/merged/ -y
done
######## DONE ####### grep -E "[A-Z0-9]{4}-10[A-Z]-[0-9]{1,3}D"
dx ls $JIE:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"*.vcf.gz" | cut -d- -f1,2,3,4 | grep -E "[A-Z0-9]{4}-10[A-Z]-01D" > blood_vardict_done
dx ls $JIE:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"*.vcf.gz" | cut -d- -f1,2,3,4 | grep -E "[A-Z0-9]{4}-10[A-Z]" > blood_mutect_done
grep -f blood_vardict_done blood_mutect_done > blood_done

for batch in {0000..0015}; do
   head -n251 dx_batch_blood.$batch.tsv > dx_batch_blood.$batch.a.tsv 
   (head -n1 dx_batch_blood.$batch.tsv; tail -n250 dx_batch_blood.$batch.tsv) > dx_batch_blood.$batch.b.tsv
   grep -v -f blood_done dx_batch_blood.$batch.a.tsv > dx_batch_blood.$batch.a2.tsv
   grep -v -f blood_done dx_batch_blood.$batch.b.tsv > dx_batch_blood.$batch.b2.tsv
done

for batch in 0008; do
   head -n251 dx_batch_blood.$batch.tsv > dx_batch_blood.$batch.a.tsv 
   (head -n1 dx_batch_blood.$batch.tsv; tail -n93 dx_batch_blood.$batch.tsv) > dx_batch_blood.$batch.b.tsv
   grep -v -f blood_done dx_batch_blood.$batch.a.tsv > dx_batch_blood.$batch.a2.tsv
   grep -v -f blood_done dx_batch_blood.$batch.b.tsv > dx_batch_blood.$batch.b2.tsv
   grep -f dx_batch_blood.$batch.a2.tsv dx_batch_blood.$batch.b2.tsv
done

######## DONE #######

dx run workflow-G7y1Fq0Jf31P502B4FF9V2Kz --batch-tsv dx_batch_blood.0005.b2.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem2_ssd1_v2_x8
sleep 1800
dx run workflow-G7y1Fq0Jf31P502B4FF9V2Kz --batch-tsv dx_batch_blood.0006.a2.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}'

251 dx_batch_blood.0005.b2.tsv
     244 dx_batch_blood.0006.a2.tsv

sleep 1800;
for batch in dx_batch_blood.0000.b2.tsv dx_batch_blood.0001.a2.tsv dx_batch_blood.0001.b2.tsv dx_batch_blood.0002.a2.tsv dx_batch_blood.0002.b2.tsv dx_batch_blood.0003.a2.tsv dx_batch_blood.0003.b2.tsv dx_batch_blood.0004.a2.tsv dx_batch_blood.0004.b2.tsv dx_batch_blood.0005.a2.tsv; do
   dx run workflow-G7y1Fq0Jf31P502B4FF9V2Kz --batch-tsv $batch -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem2_ssd1_v2_x8
done


# Tumor
# /[A-Z0-9]{4}-01[A-Z]/
grep -E "[A-Z0-9]{4}-01[A-Z]-[0-9]{2}D" tumor_done
dx ls $JIE:/CH_Exome/Workflow_Outputs/FINAL/VardictIntersectMutect/"*.vcf.gz" | cut -d- -f1,2,3,4,5 | cut -d. -f1,2 | grep -E "[A-Z0-9]{4}-01[A-Z]" | sort | uniq > tumor_vardict_done
dx ls $JIE:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/"*.vcf.gz" | cut -d- -f1,2,3,4,5 | cut -d. -f1,2 | grep -E "[A-Z0-9]{4}-01[A-Z]" | sort | uniq > tumor_mutect_done
grep -f tumor_vardict_done tumor_mutect_done > tumor_done
grep -v -f tumor_done /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/TCGA/dx_batch_tumor.0000.test.a.tsv > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/TCGA/dx_batch_tumor.0000.test.a.2.tsv

dx run workflow-G7y1Fq0Jf31P502B4FF9V2Kz --batch-tsv dx_batch_blood.0001.a2.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8
sleep 2400;
dx run workflow-G7y1Fq0Jf31P502B4FF9V2Kz --batch-tsv dx_batch_blood.0001.b2.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem2_ssd1_v2_x16
# dx ls $JIE:/CH_Exome/Workflow_Outputs/FINAL/BQSR/ | cut -d- -f1,2,3,4 | grep -E "[A-Z0-9]{4}-01[A-Z]" > tumor_bqsr_done
# grep -v -f tumor_bqsr_done /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/TCGA/dx_batch_tumor.0000.test.a.2.tsv > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/TCGA/dx_batch_tumor.0000.test.a.2.bqsr.tsv
# (head -n1 /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/TCGA/dx_batch_tumor.0000.test.a.2.tsv; grep -f tumor_bqsr_done /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/TCGA/dx_batch_tumor.0000.test.a.2.tsv) > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/TCGA/dx_batch_tumor.0000.test.a.2.finished.bqsr.tsv

# dx run workflow-G7vpZKQJf31Ppyk4BYZ9X58Q --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/TCGA/dx_batch_tumor.0000.test.a.2.bqsr.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8
# dx run workflow-G7jgVZ0Jf313xQ12J94gZ7Jp --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/TCGA/dx_batch_tumor.0000.test.a.2.finished.bqsr.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8


# dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/TCGA/dx_batch_tumor.0000.test.b.2.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8


# dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/TCGA/dx_batch_blood.0000.test.a.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem3_ssd1_v2_x8

# dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/TCGA/dx_batch_tumor.0000.test.a.tsv -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 2}}' --instance-type 0=mem2_ssd1_v2_x8


dx ls -l $JIE:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/ | grep -v TCGA > mutect.shit.list
| grep ".tbi"
tail -n+4 mutect.shit.list | while read 1 2 3 4 5 name file; do
      extension="${name##*.}"
      file=$(echo $file | awk '{gsub("\\(","",$1); print $1}' |  awk '{gsub("\\)","",$1); print $1}')
      if [[ $extension == "gz" ]]; then
         file_ext="vcf.gz"
         echo "$name   $file_ext   $file"
      else   
         file_ext="vcf.gz.tbi"
         echo "$name   $file_ext   $file"
      fi

   DX_JOB_ID=$(DD $file | grep job | awk '{print $4}')
   analysis=$(dx describe $DX_JOB_ID --json | jq -r .rootExecution)
   name_prefix=$(dx describe $analysis --json | jq -r .name | cut -d- -f2-)
   echo dx mv $file $JIE:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/$name_prefix.mutect.vep.annotated.$file_ext
   dx mv $file $JIE:/CH_Exome/Workflow_Outputs/FINAL/MutectVEP/$name_prefix.mutect.vep.annotated.$file_ext
done

 awk '{gsub("TCGA-","",$1); print}'

for file in $(head -n10 msl.txt); do 
   DX_JOB_ID=$(DD $file | grep job | awk '{print $4}')
   analysis=$(dx describe $DX_JOB_ID --json | jq -r .rootExecution)
   name_prefix=$(dx describe $analysis --json | jq -r .name)
   dx mv $file  $name_prefix.mutect.vep.annotated.vcf.gz
done


