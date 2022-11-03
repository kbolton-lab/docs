folder=40

dx build --workflow --destination $NEW_EXOME:Brian/workflows/ .
workflow=workflow-GG09ZZ8JP6z1p7f0K3B8QF5v


name=4000375
dx run $workflow -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam=file-FyK32qjJkF6115xX9VJfB62X -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index=file-FyK32qjJkF65BXqZ55P0q53P -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --name $name --instance-type 1=mem2_ssd1_v2_x8 --instance-type 3=mem2_ssd1_v2_x32

name=4072987
dx run $workflow -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam=file-G30kK58JkF68x6KG8kJ1GZYq -istage-G4Q4Pk8J6XGFjX9F3xfbg8Qv.bam_index=file-G30kzg8JkF6K14v1551XF2B4 -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --name $name --instance-type 1=mem2_ssd1_v2_x8  --instance-type 2=mem1_ssd1_v2_x8 --instance-type 3=mem2_ssd1_v2_x32


for file in $(dx ls "project-G3Yj1vjJ6XG579jbKyjXPGGY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/12/*.bqsr.HsMetrics.txt" | head -n100); do
   dx download "project-G3Yj1vjJ6XG579jbKyjXPGGY:/CH_Exome/Workflow_Outputs/FINAL/Alignment_QC/12/$file"
done

dx generate_batch_inputs -icram_in='(.*)_23193_0_0.cram$'  --path "project-G428g40Jf312vx5727kY9KXZ:/Brian/Inputs/NORMALS/" -o "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/40_WGS/40_WGS"
dx run applet-GFzyFBQJQ28JK8Py6VG913Pq --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/40_WGS/40_WGS.0000.new.tsv -iref_fasta=file-G3jKk3QJ6XG9KY8680qg9j8b -y --priority low --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 5}}' --destination project-G4qpk1jJQ285yvbXPFZKXkk8:/Brian/Workflow_Outputs/WGS/

dx ls project-G4qpk1jJQ285yvbXPFZKXkk8:/Brian/Workflow_Outputs/WGS/ | cut -d_ -f1 | sort | uniq > wgs.txt
grep -v -f wgs.txt /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/40_WGS/40_WGS.0000.tsv > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/40_WGS/40_WGS.0000.new.tsv

(printf "MEAN_TARGET_COVERAGE\tMEDIAN_TARGET_COVERAGE\tMAX_TARGET_COVERAGE\tMIN_TARGET_COVERAGE\n";
for file in $(ls *.HsMetrics.txt); do
   head -n8 $file | tail -n1 | awk '{print $34,$35,$36,$37}' OFS='\t'
done | sort -k1,1rn) | less

(printf "MEAN_TARGET_COVERAGE\tMEDIAN_TARGET_COVERAGE\tMAX_TARGET_COVERAGE\tMIN_TARGET_COVERAGE\n"; cat test.txt) | column -t 


/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/interval_lists/bed/hg38.WGS.interval_list

(grep ^@ hg38.WGS.interval_list; grep -v ^@ hg38.WGS.interval_list | awk -F'\t' ' {print $1,$2,$3,"-","."}' OFS='\t')

