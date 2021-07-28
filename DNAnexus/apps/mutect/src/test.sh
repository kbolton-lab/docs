

output_dir="out"
mkdir $output_dir
interval_list="/Users/brian/Bolton/CWL_TESTS/xgen_plus_spikein.GRCh38_chr22.interval_list"
scatter_count=10

    #bash -c "/usr/bin/java -jar /usr/picard/picard.jar IntervalListTools OUTPUT=/Users/brian/${output_dir} INPUT=${interval_list} SCATTER_COUNT=${scatter_count}"

docker run -v /Users/brian/Bolton/CWL_TESTS/:/Users/brian/Bolton/CWL_TESTS/ -v /Users/brian:/Users/brian -w /Users/brian broadinstitute/picard:2.23.6 \
    /usr/bin/perl /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/mutect/resources/usr/bin/split_interval_list_helper.pl /Users/brian/${output_dir} ${interval_list} ${scatter_count}

 docker run -v /Users/brian/Bolton/CWL_TESTS/:/Users/brian/Bolton/CWL_TESTS/ -v $HOME:$HOME -w $HOME broadinstitute/picard:2.23.6 \
    bash -c "ls ${output_dir} && echo ${output_dir}"
/Users/brian/Bolton/CWL_TESTS/xgen_plus_spikein.GRCh38_chr22.interval_list


docker run -it -v /Users/brian/Bolton/CWL_TESTS/:/Users/brian/Bolton/CWL_TESTS/ -v $HOME:$HOME -w $HOME broadinstitute/picard:2.23.6
output_dir="out"
interval_list="/Users/brian/Bolton/CWL_TESTS/xgen_plus_spikein.GRCh38_chr22.interval_list"
scatter_count=10