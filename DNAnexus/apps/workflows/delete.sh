for folder in 33; do
   for dir in FilterMutectCalls MapQ0 VardictJava; do
      dx rm $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/"*.vcf.*"
   done
done

for folder in 33; do
   for dir in MSK_Pileup; do
      dx rm $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$dir/$folder/"*.mutect.msk.pon2.vcf.*"
   done
done