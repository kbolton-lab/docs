#!/bin/bash
# Nsamples
# input:
    # column - column name to count
    # N - nsample number
    # input file
col=$1
N=$2
file=$3

column=$(COL $file $col | cut -d: -f1)

awk -F'\t' -v column="$column" '{ print $column}' $file | sort | uniq -c | sort -k1,1rn | awk -v nsamp=$N '$1 > nsamp {print}' > nsamples.remove.txt
awk '{print $1,$2'} OFS='\t' nsamples.remove.txt > nsamples.remove.txt.tmp && mv nsamples.remove.txt.tmp nsamples.remove.txt
awk -v column="$column" -F'\t' 'NR == FNR {a[$2]; next} !($column in a)' nsamples.remove.txt $file

awk -v column="$column" -F'\t' 'NR == FNR {a[$2]; next} !($column in a)' nsamples.remove.txt lung.ch_pd2.tsv > $(basename $(basename $file .txt) .tsv).nsamples.removed.tsv


awk -v column="$column" -F'\t' 'NR == FNR {a[$2]; next} !($column in a)' nsamples.remove.txt lung.hotspots.tsv > lung.hotspots.nsamples.removed.tsv

awk -v column="$column" -F'\t' 'NR == FNR {a[$2]; next} !($column in a)' nsamples.remove.txt lung.ch_pd2.tsv > lung.ch_pd2.nsamples.removed.tsv