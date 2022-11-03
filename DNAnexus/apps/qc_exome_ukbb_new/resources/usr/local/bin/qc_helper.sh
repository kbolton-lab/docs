#!/bin/bash

export bam_path=$1; shift; 
export reference_path=$1; shift; 
export reference_index_path=$1; shift;
export picard_metric_accumulation_level=$1; shift; 
echo $picard_metric_accumulation_level
export bait_intervals_path=$1; shift; 
echo $bait_intervals_path
export target_intervals_path=$1; shift; 
export per_target_coverage=$1; shift; 
export per_base_coverage=$1; shift; 
export minimum_mapping_quality=$1; shift; 
export minimum_base_quality=$1; shift; 
export alignment_summary_metrics=$1; shift; 
export hs_metrics=$1;

/usr/local/openjdk-8/bin/java -Xmx64g -jar /usr/picard/picard.jar CollectAlignmentSummaryMetrics \
    INPUT=$bam_path \
    REFERENCE_SEQUENCE=$reference_path \
    METRIC_ACCUMULATION_LEVEL=$picard_metric_accumulation_level \
    OUTPUT=$alignment_summary_metrics


/usr/local/openjdk-8/bin/java -Xmx64g -jar /usr/picard/picard.jar CollectHsMetrics \
    -I $bam_path \
    -R $reference_path \
    -METRIC_ACCUMULATION_LEVEL ALL_READS \
    -BI $bait_intervals_path \
    -TI $target_intervals_path \
    -PER_TARGET_COVERAGE $per_target_coverage \
    -PER_BASE_COVERAGE $per_base_coverage \
    -MINIMUM_MAPPING_QUALITY $minimum_mapping_quality \
    -MINIMUM_BASE_QUALITY $minimum_base_quality \
    -O $hs_metrics