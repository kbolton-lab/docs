docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin/:/usr/local/bin -v /usr/local/reference/:/usr/local/reference/ -w /home/dnanexus mgibio/mapq0-filter:v0.3.1 \
        /bin/bash /usr/bin/mapq0_vcf_filter.sh $output_name pre_mapq0.vcf.gz $tumor_bam_name $reference_name 0.15

dx download file-G4Xqyp0JYJzKfB909YjvJP0j ## pre_mapq0.vcf.gz
dx download file-G4XqzZQJj8VgfGyB8vKGQfqP ## mutect.gnomAD_AF_filter
dx download file-G4XqJg0JxZXPvP2z9YGqjk2z # $tumor_bam_name
dx download file-G4XqJqjJxZXJP9YZ8bYZYPQk

output_name=pre_gnomad/post_mapq0.vcf.gz
tumor_bam_name=1000144_23153_0_0.bqsr.bam
reference_name=GRCh38_full_analysis_set_plus_decoy_hla.fa
pre_mapq0=1000144_23153_0_0.normalized.merged.mutect.filtered.vcf.gz
bsub -oo pre_gnomad.txt -n8 -G compute-timley -g /bwileytest -q general -M 16G -R 'rusage[mem=16G]' -a 'docker(mgibio/mapq0-filter:v0.3.1)' \
/bin/bash /usr/bin/mapq0_vcf_filter.sh $output_name $pre_mapq0 $tumor_bam_name $reference_name 0.15



# arguments
outdir=$1
vcf=$2
bam=$3
mapq0perc=$4
sample_name=$5

# output_name=post_mapq0_gnomad2.vcf.gz
# post_gnomad_pre_mapq0=1000144_23153_0_0.mutect.gnomAD_AF_filter.vcf.gz
# vcf=1000051_23153_0_0.mutect.gnomAD_AF_filter.vcf.gz
# bsub -oo gnomad.txt -n8 -G compute-timley -g /bwileytest -q general -M 16G -R 'rusage[mem=16G]' -a 'docker(mgibio/mapq0-filter:v0.3.1)' \
    
# arguments
outdir=$1
vcf=$2
bam=$3
mapq0perc=$4
sample_name=$5

bsub -oo 1000051_gnomad_latest.txt -n8 -G compute-timley -g /bwileytest -q general -M 16G -R 'rusage[mem=16G]' -a 'docker(kboltonlab/mapq0-filter)' \
    /bin/bash /usr/bin/mapq0_vcf_filter.sh \
    1000051_latest \
    1000051_23153_0_0.mutect.gnomAD_AF_filter.vcf.gz \
    1000051_23153_0_0.bqsr.bam \
    UKB_4926522_233027596 \
    0.15

docker run --rm -v /Volumes/bolton/Active/projects/mocha/UKBB/exome/test:/Volumes/bolton/Active/projects/mocha/UKBB/exome/test -w /Volumes/bolton/Active/projects/mocha/UKBB/exome/test kboltonlab/mapq0-filter:v0.3.1 /bin/bash /usr/bin/test_mapq0_vcf_filter.sh \
    1000051_v0.3.1/1000051_mapq_filtered.vcf.gz \
    1000051_23153_0_0.mutect.gnomAD_AF_filter.vcf.gz \
    1000051_23153_0_0.bqsr.bam \
    GRCh38_full_analysis_set_plus_decoy_hla.fa \
    0.15

docker run --rm -v /Volumes/bolton/Active/projects/mocha/UKBB/exome/test:/Volumes/bolton/Active/projects/mocha/UKBB/exome/test -w /Volumes/bolton/Active/projects/mocha/UKBB/exome/test kboltonlab/mapq0-filter:v0.3.1 /bin/bash /usr/bin/mapq0_vcf_filter_samtools.sh \
    1000051_v0.3.1/1000051_mapq_filtered.vcf.gz \
    1000051_23153_0_0.mutect.gnomAD_AF_filter.chr17.vcf.gz \
    1000051_23153_0_0.bqsr.bam \
    GRCh38_full_analysis_set_plus_decoy_hla.fa \
    0.15

time (zgrep -v "^#" "$vcf" | grep -v "MQ0" | cut -f 1,2 | head -n100 | while read chr pos; do
    mapq0=$(samtools stats -d -@8 "$bam" $chr:$pos-$pos |  grep "reads MQ0:" | cut -f3); printf "$chr\t$pos\t$mapq0\n"
done > samtools.txt)
(printf "chr\tpos\tpysamstats\tsamtools stats\n"; paste pysamstats.txt samtools.txt | cut -f1,2,3,6 | ggrep -v -P "\t0\t0") | column -t

awk 'FNR==NR{a[FNR]=$3;next};{$NF=a[FNR]};3' samtools.txt pysamstats.txt | head

# 100: 
time (zgrep -v "^#" "$vcf" | grep -v "MQ0" | cut -f 1,2 | head -n10 | while read chr pos; do
    pysamstats --type mapq --chromosome $chr --start $pos --end $((pos+1)) "$bam"  | grep $pos | cut -f 1,2,5
done > pysamstats.txt)

chr1	1520272	216	215	14	13	50	50	60	60 # pos=1520272

# 100: 9.91s user
time (zgrep -v "^#" "$vcf" | grep -v "MQ0" | cut -f 1,2 | head -n10 | while read chr pos; do
    mapq0=$(samtools stats -@8 "$bam" $chr:$pos-$((pos+1)) |  grep "reads MQ0:" | cut -f3); printf "$chr\t$pos\t$mapq0\n"
done > samtools.txt)

chr1:1223851-1223852
# arguments
outvcf=$1
vcf=$2
bam=$3
ref_fasta=$4
mapq0perc=$5
outdir=$(dirname "$outvcf")

bsub -oo 1000051_gnomad_v0.3.1.txt -n8 -G compute-timley -g /bwileytest -q general -M 16G -R 'rusage[mem=16G]' -a 'docker(mgibio/mapq0-filter:v0.3.1)' \
    /bin/bash /usr/bin/mapq0_vcf_filter.sh \
    1000051_v0.3.1/1000051_mapq_filtered.vcf.gz \
    1000051_23153_0_0.mutect.gnomAD_AF_filter.vcf.gz \
    1000051_23153_0_0.bqsr.bam \
    GRCh38_full_analysis_set_plus_decoy_hla.fa \
    0.15






    docker run --rm -v /Volumes/bolton/Active/projects/mocha/UKBB/exome/test:/Volumes/bolton/Active/projects/mocha/UKBB/exome/test -w /Volumes/bolton/Active/projects/mocha/UKBB/exome/test kboltonlab/bst which samtools