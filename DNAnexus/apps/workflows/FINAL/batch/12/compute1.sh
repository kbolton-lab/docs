folder=12
dx download $BRIAN:/CH_Exome/Workflow_Outputs/FINAL/Final_Outputs/$folder/"*.final.tsv" --lightweight --no-progress

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

awk -F'\t' '
FNR==1 && NR!=1 { while (/^CHROM/) getline; }
    1 {print}
' OFS='\t' *.tsv > ../NSAMPS/combined.pon.germline3.tsv

function NSAMPLES {
#!/bin/bash
# Nsamples
# input:
    # column - column name to count <key or var_key>
    # N - nsample number <from MACBOOK "NSAMP N">
    # input file <combined input file>
col=$1
N=$2
file=$3

column=$(COL $file $col | cut -d: -f1)

awk -F'\t' -v column="$column" '{ print $column}' $file | sort | uniq -c | sort -k1,1rn | awk -v nsamp=$N '$1 > nsamp {print}' > nsamples.remove.txt
awk '{print $1,$2'} OFS='\t' nsamples.remove.txt > nsamples.remove.txt.tmp && mv nsamples.remove.txt.tmp nsamples.remove.txt
#awk -v column="$column" -F'\t' 'NR == FNR {a[$2]; next} !($column in a)' nsamples.remove.txt $file

awk -v column="$column" -F'\t' 'NR == FNR {a[$2]; next} !($column in a)' nsamples.remove.txt $file > $(basename $(basename $file .txt) .tsv).nsamples.removed.tsv

}


function mutect_passed () {
    awk -F "\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
        NR > 1 && ( $(f["ch_pd2"])==1 &&
                    $(f["Mutect2_PASS"])==1 &&
                    $(f["Vardict_CALLER"])==1 &&
                    $(f["alt_strand_counts_min_1_caller_only"])=="TRUE" &&
                    $(f["max.over.0.02"])=="TRUE" &&
                    ( $(f["max.under.0.35"])=="TRUE" ||
                        ($(f["max.under.0.35"])=="FALSE" && 
                    ($(f["n.loci.vep"])>=5 && $(f["n.loci.vep"])!="NA") )) &&
                    $(f["Vardict_PON_2AT2_percent"])!=1 &&
                    $(f["Mutect2_PON_2AT2_percent"])!=1 )
        ' $1 > $2
    (head -n1 $1; cat $2) > $2.tmp && mv $2.tmp $2
}
export -f mutect_passed 
mutect_passed combined.pon.germline3.nsamples.removed.tsv germline3.ch_pd2.pon.nsamples.removed.tsv
## orig
 mutect_passed mutect_passed combined.pon.germline2.tsv germline2.ch_pd2.pon.nsamples.removed.NOT.stream_gene_variant.tsv # 1657 lines


function hotspots_called () {
    awk -F "\t" 'NR==1 {for (i=1; i<=NF; i++) f[$i] = i}
        NR > 1 && ( $(f["called"])=="TRUE" &&
                    $(f["Mutect2_PASS"])!=1 &&
                    ( $(f["heme_cosmic_count"])>=10 || $(f["myeloid_cosmic_count"])>=5 || 
                    ($(f["n.loci.vep"])>=5 && $(f["n.loci.vep"])!="NA") || $(f["n_panmyeloid"])>=5 ) &&
                    $(f["Vardict_PON_2AT2_percent"])!=1 &&
                    $(f["Mutect2_PON_2AT2_percent"])!=1 &&
                    $(f["alt_strand_counts_min_1_caller_only"])=="TRUE" )
        ' $1 > $2
    (head -n1 $1; cat $2) > $2.tmp && mv $2.tmp $2
}
export -f hotspots_called
hotspots_called combined.pon.germline2.nsamples.removed.tsv germline2.hotspots.pon.nsamples.removed.tsv