function TOTAL {
    for subfolder in {12..13};
        # K
        dx ls UKBB_Exome_2021:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/$subfolder | cut -d_ -f1 | sort | uniq >>  /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/brian_project
        cat /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/brian_project | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/brian_project.tmp && mv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/brian_project.tmp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/brian_project

        dx ls Kelly_UKBB:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/$subfolder | cut -d_ -f1 | sort | uniq >>  /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/kelly_project
        cat /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/kelly_project | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/kelly_project.tmp && mv /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/kelly_project.tmp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/kelly_project

        cat /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/brian_project /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/kelly_project | sort | uniq > /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/TOTAL/total.txt
}