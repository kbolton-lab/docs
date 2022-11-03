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

    set -ex -o pipefail
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
    echo "${dockerimage_picard_path}"
    
    my_python=$(which python3)
    echo $my_python
    multiqc -h
    # ln -s $my_python /usr/bin/python3
    
    
    dx-download-all-inputs --parallel
    mv $reference_index_path $(dirname $reference_path)
    ls $(dirname $reference_path)
    per_target_coverage=${bam_prefix}.PerTargetCoverage.txt
    per_base_coverage=${bam_prefix}.PerBaseCoverage.txt
    alignment_summary_metrics=${bam_prefix}.AlignmentSummaryMetrics.txt
    hs_metrics=${bam_prefix}.HsMetrics.txt

    docker load -i ${dockerimage_picard_path}

    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin:/usr/local/bin -w /home/dnanexus broadinstitute/picard:2.23.6 /usr/local/bin/qc_helper.sh $bam_path $reference_path $reference_index_path $picard_metric_accumulation_level $bait_intervals_path $target_intervals_path $per_target_coverage $per_base_coverage $minimum_mapping_quality $minimum_base_quality $alignment_summary_metrics $hs_metrics

    multiqc .
    ls multiqc_data
    mv multiqc_data/multiqc_data.json "${bam_prefix}.multiqc_data.json"
    mv multiqc_data/multiqc_general_stats.txt "${bam_prefix}.multiqc_general_stats.txt"

    alignment_summary_metrics=$(dx upload $alignment_summary_metrics --brief)
    hs_metrics=$(dx upload $hs_metrics --brief)
    per_target_coverage=$(dx upload $per_target_coverage --brief)
    multiqc_data=$(dx upload "${bam_prefix}.multiqc_data.json" --brief)
    multiqc_general_stats=$(dx upload "${bam_prefix}.multiqc_general_stats.txt" --brief)
    #per_base_coverage=$(dx upload $per_base_coverage --brief)


    dx-jobutil-add-output alignment_summary_metrics "$alignment_summary_metrics" --class=file
    dx-jobutil-add-output hs_metrics "$hs_metrics" --class=file
    dx-jobutil-add-output multiqc_data "$multiqc_data" --class=file
    dx-jobutil-add-output multiqc_general_stats "$multiqc_general_stats" --class=file
    #dx-jobutil-add-output per_base_coverage "$per_base_coverage" --class=file
}