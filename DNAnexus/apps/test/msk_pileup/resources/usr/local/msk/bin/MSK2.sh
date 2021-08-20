reference=$1
bams=$2
filename=$3
msk_out=$4

echo "/opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta $reference ${bams} --vcf $filename --output $msk_out --thread 8;"
echo "bgzip -f $msk_out && tabix -f $msk_out.gz"