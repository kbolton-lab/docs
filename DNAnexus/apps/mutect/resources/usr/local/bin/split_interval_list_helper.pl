use File::Copy;

die "wrong number of inputs" unless scalar(@ARGV) == 3;
my ($output_dir, $interval_list, $scatter_count) = @ARGV;

my $i = 1;

if ($scatter_count == 1) {
    File::Copy::copy($interval_list,qq{$i.interval_list});
} else {

    my $retval = system('/usr/bin/java', '-jar', '/usr/picard/picard.jar', 'IntervalListTools', 'OUTPUT='.$output_dir, 'INPUT='.$interval_list, 'SCATTER_COUNT='. $scatter_count);
    exit $retval if $retval != 0;

    for (glob("$output_dir/*/scattered.interval_list")) {
        #create unique names and relocate all the scattered intervals to a single directory
        File::Copy::move($_, qq{$output_dir/$i.interval_list});
        $i++
    }
}