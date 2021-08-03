docker run --rm -v /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/mutect/test:/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/mutect/test -v /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/mutect/resources/usr/bin:/usr/local/bin -v /Users/brian:/Users/brian -w /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/mutect/test broadinstitute/picard:2.23.6 \
        /usr/bin/perl /usr/local/bin/split_interval_list_helper.pl /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/mutect/test/${output_dir} ${input_to_scatter} ${scatter_count}

scatter() {
    echo "Value of input_to_scatter: '${input_to_scatter}'"

    # Fill in code here to do whatever is necessary to scatter the
    # input.
    docker load -i /picard.tar.gz

    output_dir="out"
    mkdir $output_dir
    scatter_count=10

    docker run --rm -v /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/mutect/test:/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/mutect/test -v /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/mutect/resources/usr/bin:/usr/local/bin -v /Users/brian/Bolton/CWL_TESTS/:/Users/brian/Bolton/CWL_TESTS/ -w /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/mutect/test broadinstitute/picard:2.23.6 \
        /usr/bin/perl /usr/local/bin/split_interval_list_helper.pl /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/mutect/test/${output_dir} ${input_to_scatter} ${scatter_count}

    #declare -a scattered_input=(placeholder1 placeholder2)
    scattered_input=( $(ls ${output_dir}/*.interval_list) )

    for piece in "${scattered_input[@]}"
    do
        dx-jobutil-add-output array_of_scattered_input "$piece" --array
    done
}

mkdir /tmp/mutect
cd /tmp/mutect
cp /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/mutect/test/out/*.interval_list .
ls -1sh /tmp/mutect/*.interval_list