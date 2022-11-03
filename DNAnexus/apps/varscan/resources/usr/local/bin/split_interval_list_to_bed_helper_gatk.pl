#!/usr/bin/perl

use strict;
use warnings;

use feature qw(say);
use Data::Dumper;
#my $retval = system('gatk', 'IntervalListTools', @ARGV);
#exit $retval if $retval != 0;

my $i = 1;
# my $XMLDIR = ".";
#https://stackoverflow.com/questions/25206041/sorting-file-names-by-numeric-value
my @files = sort {(split /\./, $a)[0] <=> (split /\./, $b)[0]} glob("*.interval_list");
print Dumper(@files2);


for my $interval_list (@files) {
    my $bed = $i.'.interval.bed';
    open(my $in_fh, $interval_list) or die "fail to open $interval_list for read"; 
    open(my $out_fh, ">$bed") or die "fail to write to $bed";
    while (<$in_fh>) {
        next if /^@/;
        my ($chr, $start, $stop) = split /\t/, $_;
        $out_fh->say(join("\t", $chr, $start-1, $stop));
    }
    close $in_fh;
    close $out_fh;
    $i++
}