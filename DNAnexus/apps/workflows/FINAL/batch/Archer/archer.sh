i=1
REF=/Volumes/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa
for bam in $(ls /Volumes/bolton/Active/projects/Dilution/hg38/Downsample/272341_1*.bam) /Volumes/bolton/Active/projects/Dilution/hg38/Downsample/736399_T_S221.bam; do
    
    DIR=/Volumes/bolton/Active/projects/Dilution/hg38/Downsample
    name=$(basename $bam .bam)
    samtools view $bam -L $HG38/bed/archer_panel.hg38.sorted.bed -o $DIR/brian/${name}_bedtargets.bam
    samtools index $DIR/brian/${name}_bedtargets.bam

    java -jar /Users/brian/tools/jvarkit/dist/sortsamrefname.jar \
        --samoutputformat BAM \
        $DIR/brian/${name}_bedtargets.bam |\
    java -jar /Users/brian/tools/jvarkit/dist/biostar154220.jar \
        -d 55 \
        --regions $HG38/bed/archer_panel.hg38.sorted.bed \
        --samoutputformat BAM |\
    samtools sort -o $DIR/brian/${name}_ukbbPipe_sample_${i}.bam
    samtools index $DIR/brian/${name}_ukbbPipe_sample_${i}.bam
    
    dx upload $DIR/brian/${name}_ukbbPipe_sample_${i}.bam* --destination  $KELLY:/CH_Exome/Inputs/Archer/jvarkit/bams_${i}/
done




workflow=workflow-G5kgvY0JQ28BYfXQ8Q0z8Q61
function dilution_jvarkit {
    sample=$1
    normal=$2
    workflow=$3
    wf_folder=/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/Archer

    for i in {1..1}; do
        for j in 10 100 1000 5000; do
            tumor_bam=$(dx find data --class file --path $KELLY:/CH_Exome/Inputs/Archer/jvarkit/bams_${i}/ --json --name ${sample}_1_${j}_ukbbPipe_sample_${i}.bam | jq -c '.[].describe.id')
            tumor_bai=$(dx find data --class file --path $KELLY:/CH_Exome/Inputs/Archer/jvarkit/bams_${i}/ --json --name ${sample}_1_${j}_ukbbPipe_sample_${i}.bam.bai | jq -c '.[].describe.id')

            normal_bam=$(dx find data --class file --path $KELLY:/CH_Exome/Inputs/Archer/jvarkit/bams_${i}/ --json --name ${normal}_ukbbPipe_sample_${i}.bqsr.bam | jq -c '.[].describe.id')
            normal_bai=$(dx find data --class file --path $KELLY:/CH_Exome/Inputs/Archer/jvarkit/bams_${i}/ --json --name ${normal}_ukbbPipe_sample_${i}.bqsr.bam.bai | jq -c '.[].describe.id')
            
            /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/archer_in.py ${sample}_1_${j}.json $tumor_bam $tumor_bai $normal_bam $normal_bai $wf_folder
            
            echo dx run $workflow -f $wf_folder/${sample}_1_${j}.json -y --priority low --project $KELLY --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 3}}'
        done
    done
}
export -f dilution_jvarkit
dilution_jvarkit 272341 736399_T_S221 $workflow

workflow=workflow-G5kgvY0JQ28BYfXQ8Q0z8Q61
function dilution_brian {
    sample=$1
    normal=$2
    workflow=$3
    wf_folder=/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/Archer

    for i in {1..1}; do
        for j in 10; do
            tumor_bam=$(dx find data --class file --path $KELLY:/CH_Exome/Inputs/Archer/ --json --name ${sample}_1_${j}.test.sorted.bam | jq -c '.[].describe.id')
            tumor_bai=$(dx find data --class file --path $KELLY:/CH_Exome/Inputs/Archer/ --json --name ${sample}_1_${j}.test.sorted.bam.bai | jq -c '.[].describe.id')

            normal_bam=$(dx find data --class file --path $KELLY:/CH_Exome/Inputs/Archer/ --json --name ${normal}.test.sorted.bqsr.bam | jq -c '.[].describe.id')
            normal_bai=$(dx find data --class file --path $KELLY:/CH_Exome/Inputs/Archer/ --json --name ${normal}.test.sorted.bqsr.bam.bai | jq -c '.[].describe.id')
            
            /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/archer_in.py ${sample}_1_${j}.json $tumor_bam $tumor_bai $normal_bam $normal_bai $wf_folder
            
            echo dx run $workflow -f $wf_folder/${sample}_1_${j}.json -y --priority low --project $KELLY --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 3}}'
        done
    done
}
export -f dilution_brian
dilution_brian 272341 736399_T_S221 $workflow

dx run workflow-G5kgvY0JQ28BYfXQ8Q0z8Q61 -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/Archer/272341_1_10.json -y --priority low --project project-G4qpk1jJQ285yvbXPFZKXkk8 --extra-args '{"executionPolicy": {"restartOn": {"UnresponsiveWorker": 2}, "maxSpotTries": 3}}'


awk -F'\t' '
FNR==1 && NR!=1 { while (/^GENOMIC/) getline; }
    1 {print}
' OFS='\t' CosmicMutantExport.final.high.myeloid.chr*.tsv > CosmicMutantExport.high.myeloid.tsv

awk -F'\t' '
FNR==1 && NR!=1 { while (/^GENOMIC/) getline; }
    1 {print}
' OFS='\t' CosmicMutantExport.final.high.heme.chr*.tsv > CosmicMutantExport.high.heme.tsv

DIR=/storage1/fs1/bolton/Active/projects/Dilution/hg38/Downsample
 # depth=$(samtools depth -b $HG38/bed/archer_panel.hg38.sorted.bed $DIR/brian/${name}_bedtargets.bam | AVG 3)
samtools depth -b ASXL.bed $DIR/brian/${name}_bedtargets.bam | AVG 3
    # RATIO_OF_DOWNSAMPLING=$(echo "scale=4; 56/$depth" | bc)

touch test.ASXL1.exon1.bam
cat ASXL_single.bed | while read chr start end; do
    echo samtools view -r $chr:$start-$end -@16 -b $DIR/brian/${name}_bedtargets.bam
    samtools view -r $chr:$start-$end -@16 -b $DIR/brian/${name}_bedtargets.bam | wc -l
done

bam=272341_1_10.bam
name=$(basename $bam .bam)
DIR=/storage1/fs1/bolton/Active/projects/Dilution/hg38/Downsample
rm test.asx1l.exon1.sam
touch test.asx1l.exon1.sam
samtools view -H $DIR/brian/${name}_bedtargets.bam -b > header.txt
cat ASXL_single.bed | while read chr start end; do
    #count=$(samtools view -@16 -c $DIR/brian/${name}_bedtargets.bam $chr:$start-$end)
    samtools view -@16 $DIR/brian/${name}_bedtargets.bam $chr:$start-$end | shuf | head -n56 >> test.asx1l.exon1.sam
    #depth=$(samtools depth -r $chr:$start-$end $DIR/brian/${name}_bedtargets.bam -J | awk '{print $3}')
    # if [[ $count -eq $depth ]]; then
    #     echo yes
    #     $chr:$start-$end
    #     echo $depth
    #     echo $count
    # else 
    #     echo no
    #     $chr:$start-$end
    #     echo $depth
    #     echo $count
    # fi
done
(samtools view header.txt -h; sort test.asx1l.exon1.sam | uniq) | samtools depth
samtools depth  -J | awk '{print $3}'

bam=736399_T_S221.bam
bam=272341_1_10.bam
name=$(basename $bam .bam)
samtools bedcov ASXL_single.bed ${name}_bedtargets.bam | awk -F'\t' '$5=$4/($3-$2) {print}' OFS='\t' | awk -F'\t' '$6=56/$5 {print}' OFS='\t' > archer.exon.depth.avg.$name.bed

samtools view -H $DIR/brian/${name}_bedtargets.bam > $name.test.bam
cat archer.exon.depth.avg.$name.bed | while read chr start end total avg perc; do
    count=$(samtools view -@16 $DIR/brian/${name}_bedtargets.bam $chr:$start-$end -c) 
    lines=$(printf "%.0f" $(echo $perc*$count | bc))
    printf "$count\t$lines\n" 
    samtools view -@16 $DIR/brian/${name}_bedtargets.bam $chr:$start-$end | shuf | head -n $lines >> $name.test.bam
    sort $name.test.bam | uniq > $name.test.bam.tmp && mv $name.test.bam.tmp $name.test.bam
done
samtools sort $name.test.bam -o $name.test.sorted.bam && samtools index $name.test.sorted.bam
samtools bedcov ASXL_single.bed $name.test.sorted.bam | awk -F'\t' '$5=$4/($3-$2) {print}' OFS='\t'
samtools bedcov $HG38/bed/archer_panel.hg38.sorted.bed $name.test.sorted.bam | awk -F'\t' '$6=$5/($3-$2) {print}' OFS='\t'




    mapq0=$(samtools stats -d -@8 "$bam" $chr:$pos-$pos | grep "reads MQ0:" | cut -f3); printf "$chr\t$pos\t$mapq0\n" >> mapq0counts

    /gatk/gatk DownsampleSam \
        -I $DIR/brian/${name}_bedtargets.bam \
        -O $DIR/brian/${name}_ukbbPipe_sample_${i}_gatk.bam \
        -S HighAccuracy \
        -P $RATIO_OF_DOWNSAMPLING
    
    # -

    mkdir /tmp/spark
    gatk --java-options '-Xmx32G -XX:+UseParallelGC -XX:ParallelGCThreads=8' BQSRPipelineSpark \
        -R $GMSROOT/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa \
        -I /Volumes/bolton/Active/projects/Dilution/hg38/Downsample/brian/218281_1_10_ukbbPipe_sample_1.bam \
        -O /Volumes/bolton/Active/projects/Dilution/hg38/Downsample/brian/218281_1_10_ukbbPipe_sample_1.bqsr.bam \
        --known-sites $GMSROOT/gc2560/core/build_merged_alignments/detect-variants--linus2112.gsc.wustl.edu-jwalker-20211-26b393cc7ab04120ac68cc2cbd4a15df/indels.hq.vcf.gz \
        --known-sites $GMSROOT/gc2560/core/build_merged_alignments/detect-variants--linus2112.gsc.wustl.edu-jwalker-19443-e48c595a620a432c93e8dd29e4af64f2/snvs.hq.vcf.gz
        --verbosity ERROR -- --spark-runner LOCAL --spark-master 'local[8]' --conf spark.local.dir=/tmp/spark

    gatk --java-options '-Xmx32G -XX:+UseParallelGC -XX:ParallelGCThreads=8' BQSRPipelineSpark \
        -R $GMSROOT/gc2560/core/model_data/2887491634/build21f22873ebe0486c8e6f69c15435aa96/all_sequences.fa \
        -I /Volumes/bolton/Active/projects/Dilution/hg38/Downsample/218281_1_10.bam \
        -O /Volumes/bolton/Active/projects/Dilution/hg38/Downsample/brian/218281_1_10.bqsr.bam \
        --known-sites $GMSROOT/gc2560/core/build_merged_alignments/detect-variants--linus2112.gsc.wustl.edu-jwalker-20211-26b393cc7ab04120ac68cc2cbd4a15df/indels.hq.vcf.gz \
        --known-sites $GMSROOT/gc2560/core/build_merged_alignments/detect-variants--linus2112.gsc.wustl.edu-jwalker-19443-e48c595a620a432c93e8dd29e4af64f2/snvs.hq.vcf.gz
        --verbosity ERROR -- --spark-runner LOCAL --spark-master 'local[8]' --conf spark.local.dir=/tmp/spark


    /Volumes/bolton/Active/projects/Dilution/hg38/Downsample/brian/218281_1_10_ukbbPipe_sample_1.bam