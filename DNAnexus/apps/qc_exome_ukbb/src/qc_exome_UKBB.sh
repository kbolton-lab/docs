#!/bin/bash
# qc_exome_UKBB 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://documentation.dnanexus.com/developer for tutorials on how
# to modify this file.

main() {

    echo "Value of bam: '$bam'"
    echo "Value of reference: '$reference'"
    echo "Value of reference_index: '$reference_index'"
    echo "Value of bait_intervals: '$bait_intervals'"
    echo "Value of target_intervals: '$target_intervals'"
    echo "Value of per_base_intervals: '$per_base_intervals'"
    echo "Value of omni_vcf: '$omni_vcf'"
    echo "Value of per_target_intervals: '$per_target_intervals'"
    echo "Value of picard_metric_accumulation_level: '$picard_metric_accumulation_level'"
    echo "Value of minimum_mapping_quality: '$minimum_mapping_quality'"
    echo "Value of minimum_base_quality: '$minimum_base_quality'"
    echo "${bam_prefix}"
    
    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".

    # dx download "$bam" -o bam

    # dx download "$reference" -o reference

    # dx download "$bait_intervals" -o bait_intervals

    # dx download "$target_intervals" -o target_intervals

    # dx download "$per_base_intervals" -o per_base_intervals

    # dx download "$omni_vcf" -o omni_vcf

    # dx download "$per_target_intervals" -o per_target_intervals

    dx-download-all-inputs --parallel
    mv $reference_index_path $(dirname $reference_path)
    ls $(dirname $reference_path)
    per_target_coverage="${bam_prefix}-PerTargetCoverage.txt"
    per_base_coverage="${bam_prefix}-PerBaseCoverage.txt"
    alignment_summary_metrics=${bam_prefix}.AlignmentSummaryMetrics.txt
    hs_metrics=${bam_prefix}.HsMetrics.txt

    docker load -i /picard.tar.gz

    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/bin:/usr/bin -w /home/dnanexus broadinstitute/picard:2.23.6 /usr/bin/qc_helper.sh $bam_path $reference_path $reference_index_path $picard_metric_accumulation_level $bait_intervals_path $target_intervals_path $per_target_coverage $per_base_coverage $minimum_mapping_quality $minimum_base_quality $alignment_summary_metrics $hs_metrics

    # Fill in your application code here.
    #
    # To report any recognized errors in the correct format in
    # $HOME/job_error.json and exit this script, you can use the
    # dx-jobutil-report-error utility as follows:
    #
    #   dx-jobutil-report-error "My error message"
    #
    # Note however that this entire bash script is executed with -e
    # when running in the cloud, so any line which returns a nonzero
    # exit code will prematurely exit the script; if no error was
    # reported in the job_error.json file, then the failure reason
    # will be AppInternalError with a generic error message.

    # The following line(s) use the dx command-line tool to upload your file
    # outputs after you have created them on the local file system.  It assumes
    # that you have used the output field name for the filename for each output,
    # but you can change that behavior to suit your needs.  Run "dx upload -h"
    # to see more options to set metadata.

    alignment_summary_metrics=$(dx upload $alignment_summary_metrics --brief)
    hs_metrics=$(dx upload $hs_metrics --brief)
    per_target_coverage=$(dx upload $per_target_coverage --brief)
    per_base_coverage=$(dx upload $per_base_coverage --brief)

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    dx-jobutil-add-output alignment_summary_metrics "$alignment_summary_metrics" --class=file
    dx-jobutil-add-output hs_metrics "$hs_metrics" --class=file
    dx-jobutil-add-output per_target_coverage "$per_target_coverage" --class=file
    dx-jobutil-add-output per_base_coverage "$per_base_coverage" --class=file
}
