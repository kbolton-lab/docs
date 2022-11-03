folder=10
dx download $KELLY:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/$folder/"*.final.tsv" --lightweight --no-progress

for file in $(ls *.final.tsv); do
   head -n1 $file | sha256sum
done

function COMBINE {
#!/bin/bash
# combines files with same suffix and have headers
# 1=header start
# 2=suffix
awk -F'\t' '
FNR==1 && NR!=1 { while (/^CHROM/) getline; }
    1 {print}
' OFS='\t' *.final.tsv > combined.UKBB.10.tsv
}

function NSAMPLES {
   #!/usr/bin/bash
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
#awk -v column="$column" -F'\t' 'NR == FNR {a[$2]; next} !($column in a)' nsamples.remove.txt $file

awk -v column="$column" -F'\t' 'NR == FNR {a[$2]; next} !($column in a)' nsamples.remove.txt $file > $(basename $(basename $file .txt) .tsv).nsamples.removed.tsv
}
