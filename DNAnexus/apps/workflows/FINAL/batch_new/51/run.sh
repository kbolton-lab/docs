dx run applet-G6F0v98JQ28F5pFXP2qjyqp8 -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/mutect_single/mutect_test_inputs.json --instance-type mem2_ssd1_v2_x4 --destination $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Mutect_Rerun -y

folder=51 
dx generate_batch_inputs -istage-G4KkV88J6XG340022vxq916x.tumor_bam='(.*)_23153_0_0.cram$' -istage-G4KkV88J6XG340022vxq916x.tumor_bai='(.*)_23153_0_0.cram.crai' --path "$BRIAN:Bulk/Exome sequences/Exome OQFE CRAM files/$folder" -o "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_rerun"

dx build --workflow . --destination project-G3Yj1vjJ6XG579jbKyjXPGGY:/CH_Exome/test/ 
workflow=workflow-GGb1PF8J6XGB3bZyJ966g4vq

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_rerun.0000.tsv -y --priority low --destination $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Mutect_Rerun --instance-type 0=mem2_ssd1_v2_x4

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_rerun.0001.tsv -y --priority low --destination $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Mutect_Rerun --instance-type 0=mem2_ssd1_v2_x4

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_rerun.0002.tsv -y --priority low --destination $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Mutect_Rerun --instance-type 0=mem2_ssd1_v2_x4

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_rerun.0003.tsv -y --priority low --destination $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Mutect_Rerun --instance-type 0=mem2_ssd1_v2_x4
sleep 7200
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_rerun.0004.tsv -y --priority low --destination $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Mutect_Rerun --instance-type 0=mem2_ssd1_v2_x4

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_rerun.0005.tsv -y --priority low --destination $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Mutect_Rerun --instance-type 0=mem2_ssd1_v2_x4

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_rerun.0006.tsv -y --priority low --destination $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Mutect_Rerun --instance-type 0=mem2_ssd1_v2_x4

dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_rerun.0007.tsv -y --priority low --destination $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Mutect_Rerun --instance-type 0=mem2_ssd1_v2_x4
dx run $workflow --batch-tsv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_rerun.0008.tsv -y --priority low --destination $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Mutect_Rerun --instance-type 0=mem2_ssd1_v2_x4



folder=37
dx generate_batch_inputs -istage-G4KkV88J6XG340022vxq916x.tumor_bam='(.*)_23153_0_0.cram$' -istage-G4KkV88J6XG340022vxq916x.tumor_bai='(.*)_23153_0_0.cram.crai' --path "$BRIAN:Bulk/Exome sequences/Exome OQFE CRAM files/$folder" -o "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_rerun"
dx ls "project-G3Yj1vjJ6XG579jbKyjXPGGY:/CH_Exome/Workflow_Outputs/FINAL/Mutect_Rerun/$folder/*.vcf.gz" | cut -d_ -f1 > rerun.ran

for batch in {0000..0007}; do;
    grep -v -f rerun.ran /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_rerun.$batch.tsv > tmp && mv tmp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_rerun.$batch.tsv
done
wc -l /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch_new/$folder/dx_rerun*