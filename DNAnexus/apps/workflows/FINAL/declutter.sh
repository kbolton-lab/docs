for directory in MapQ0 MSK_Pileup VardictJava; do
   for folder in 10; do 
      dx mkdir $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$directory/$folder/declutter/
   done
done

for directory in FilterMutectCalls MapQ0 MSK_Pileup VardictJava; do
   for folder in {10..60} Germline Lung; do 
      dx mkdir $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/$directory/$folder/declutter/
   done
done

project-G4qpk1jJQ285yvbXPFZKXkk8:/CH_Exome/Workflow_Outputs/FINAL/FilterMutectCalls/10/



Delete all HaplotypeCaller stuff?

dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/FilterMutectCalls/declutter/
dx ls $KELLY:/CH_Exome/Workflow_Outputs/FINAL/$directory/$folder/declutter/

$KELLY
10 12 14 16 18 20 22 24                                
