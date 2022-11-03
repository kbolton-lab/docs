sample="$(sed -n "${LSB_JOBINDEX}p" /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/${folder}/vep_run_again.txt)"
echo $sample
vcf_path=${sample}_23153_0_0.normalized.mutect.filtered.vcf.gz
# for vcf_path in $(ls *_23153_0_0.normalized.mutect.filtered.vcf.gz | grep -v 1001489 | tail -n+2 | grep -v 1000043); do
base=$(basename $vcf_path .vcf.gz)
# bsub -oo logs/$base.vep3.log -G compute-timley -g /bolton/bwileytest -q general -M 32G -R 'rusage[mem=32G]' -a 'docker(kboltonlab/ic_vep)' bash -c \
/opt/vep/src/ensembl-vep/vep \
--format vcf \
-i $vcf_path \
--transcript_version \
--offline \
--cache \
--symbol \
--vcf \
-o $base.vep.vcf \
--fasta /scratch1/fs1/bolton/brian/GRCh38_full_analysis_set_plus_decoy_hla.fa \
--dir_plugins /storage1/fs1/bolton/Active/projects/mocha/UKBB/ukbb_calls/pvcf/plugins \
--dir_cache /storage1/fs1/bolton/Active/projects/mocha/UKBB/ukbb_calls/pvcf/vep_zip \
--synonyms /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build50f99e75d14340ffb5b7d21b03887637/chromAlias.ensembl.txt \
--plugin Frameshift \
--plugin Wildtype \
--assembly GRCh38 \
--cache_version 104 \
--species homo_sapiens \
--sift p --polyphen p --pick --pick_order canonical,rank,mane,ccds,appris,tsl,biotype,length --everything 1 --merged --check_existing --buffer_size 1000 --af_gnomad \
--custom /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/genome-db-ensembl-gnomad/2dd4b53431674786b760adad60a29273/fixed_b38_exome.vcf.gz,gnomADe,vcf,exact,1,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS \
--custom /storage1/fs1/bolton/Active/data/hg38/vcf/gnomad.genomes.r3.0.sites.full_trimmed_info.af_only.vcf.gz,gnomADg,vcf,exact,1,AF,AF_ami,AF_oth,AF_afr,AF_sas,AF_asj,AF_fin,AF_amr,AF_nfe,AF_eas \
--custom /storage1/fs1/bga/Active/gmsroot/gc2560/core/custom_clinvar_vcf/v20181028/custom.vcf.gz,clinvar,vcf,exact,1,CLINSIGN,PHENOTYPE,SCORE,RCVACC,TESTEDINGTR,PHENOTYPELIST,NUMSUBMIT,GUIDELINES \
--force_overwrite  && bgzip $base.vep.vcf && tabix $base.vep.vcf.gz
# done

 PCA.r1.TCGAbarcode.merge.tnSwapCorrected.10389.vcf.gz
    
wget ftp://ftptrace.ncbi.nlm.nih.gov/1000genomes/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

ls *.vcf.gz | grep -v vep | cut -d_ -f1 > vep.run
ls *.vep.vcf.gz | cut -d_ -f1 > vep.ran
grep -v -f vep.ran vep.run > vep_run_again.txt

bsub -n16 -o bcftools.log -G compute-timley -g /bwileytest2 -q general -M 64G -R 'select[mem>64G] rusage[mem=64G]' -a 'docker(kboltonlab/bst)' /bin/bash -c "bcftools merge -m all *vep.vcf.gz -Oz -o merged/merged.vcf.gz --threads 64 && tabix merged/merged.vcf.gz"

bsub -o download.log -G compute-timley -g /bwileytest2 -q general -M 96G -R 'select[mem>96G] span[hosts=1] rusage[mem=96G]' -a 'docker(kboltonlab/dnanexus:1.0)' $EXOME/results/download.rerun.sh 28

for folder in {47..59}; do
echo $folder;
cd /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/$folder
wc -l vep_run_again.txt
ls *.vcf.gz | grep -v vep | cut -d_ -f1 > vep.run
ls *.vep.vcf.gz | cut -d_ -f1 > vep.ran
grep -v -f vep.ran vep.run > vep_run_again.txt
wc -l vep_run_again.txt
echo
done
bsub -o logs/vep_${folder}_%J_%I.out -G compute-timley -g /bwileytest2 -q general -M 4G -R 'select[mem>4G] span[hosts=1] rusage[mem=4G]' -a 'docker(kboltonlab/ic_vep)' -J "ann_${folder}[1-2]" /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/vep.sh ${folder}
# echo
# done
# bsub -o logs/vep_${folder}_%J_%I.out -G compute-timley -g /bwileytest2 -q general -M 4G -R 'select[mem>4G] span[hosts=1] rusage[mem=4G]' -a 'docker(kboltonlab/ic_vep)' -J "ann_${folder}[1001-2000]" /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/vep.sh ${folder}
# bsub -o logs/vep_${folder}_%J_%I.out -G compute-timley -g /bwileytest2 -q general -M 4G -R 'select[mem>4G] span[hosts=1] rusage[mem=4G]' -a 'docker(kboltonlab/ic_vep)' -J "ann_${folder}[2001-3000]" /storage1/fs1/bolton/Active/projects/mocha/UKBB/exfoldome/results/RERUN/vep.sh ${folder}
# bsub -o logs/vep_${folder}_%J_%I.out -G compute-timley -g /bwileytest2 -q general -M 4G -R 'select[mem>4G] span[hosts=1] rusage[mem=4G]' -a 'docker(kboltonlab/ic_vep)' -J "ann_${folder}[3001-3950]" /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/vep.sh ${folder}

for dir in {20..40}; do
    echo $dir
    cd /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/$dir
    ls merged/final/*.vcf.gz | head -n1
    echo
done



# bsub -oo $dir.PON.log -n16 -G compute-bolton -g /bwileytest -q general -M 32G -R 'select[mem>64G] rusage[mem=64]' -a 'docker(kboltonlab/msk_getbasecounts:3.1)' /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/PON.sh

bsub -oo $dir.R.log -n4 -G compute-bolton -g /bwileytest -q general -M 6G -R 'select[mem>6G] rusage[mem=6G]' -a 'docker(kboltonlab/r_bcftools)' /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/runR.sh
bsub -oo $dir.R.log -n4 -G compute-bolton -g /bwileytest -q general -M 6G -R 'select[mem>6G] rusage[mem=6G]' -a 'docker(kboltonlab/r_bcftools)' ./runR.sh

for dir in $(ls /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/*/merged/merged.norm.vcf | cut -d/ -f12 | grep -v 10 | grep -v 11); do 
    echo "################################################################################"
    echo "##################################  $dir  ###################################"
    echo "################################################################################"
    echo $dir >> /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/annotated.txt
    cd /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/$dir/merged/
    bcftools view -s $(bcftools query -l merged.norm.vcf | head -n1) merged.norm.vcf | bcftools annotate --threads 4 -a /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/normal2.txt.gz -h pon2.header -c CHROM,POS,REF,ALT,PON_2AT2_percent -Oz -o single.merged.norm.vcf.gz
    /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/annotate_pd_scratch.sh $dir
    echo
    echo
    echo "################################################################################"
    echo "##################################  $dir DONE  ##############################"
    echo "################################################################################"
    echo
    echo
done

dir=10
for sample in $(head -n1 /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/$dir/merged/final.run); do
    echo "################################################################################"
    echo "##################################  $sample  ###################################"
    echo "################################################################################"
    LSB_JOBINDEX=1
    /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/annotate_pd_scratch_SAMPLE.sh $dir $sample
    echo
    echo
    echo "################################################################################"
    echo "##################################  $sample DONE  ##############################"
    echo "################################################################################"
    echo
    echo
done

ls outputs/*.final.tsv | cut -d/ -f2 | cut -d. -f1 > final.ran
grep -v -f final.ran ../vep.run > final.run
wc -l final.run
folder=14
for sample in $(cat final.run); do
    /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/annotate_pd_scratch_SAMPLE.sh 

folder=28
bsub -o logs/final_${folder}_%J_%I.out -G compute-timley -g /bwileytest2 -q general -M 2G -R 'select[mem>2G] span[hosts=1] rusage[mem=2G]' -a 'docker(kboltonlab/annotate_wes_ch:3.1)' -J "ann_${folder}[1-2]" /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/annotate_pd_scratch_SAMPLE.sh ${folder}
bsub -o logs/final_${folder}_%J_%I.out -G compute-timley -g /bwileytest2 -q general -M 2G -R 'select[mem>2G] span[hosts=1] rusage[mem=2G]' -a 'docker(kboltonlab/annotate_wes_ch:3.1)' -J "ann_${folder}[1001-2000]" /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/annotate_pd_scratch_SAMPLE.sh ${folder}
bsub -o logs/final_${folder}_%J_%I.out -G compute-timley -g /bwileytest2 -q general -M 2G -R 'select[mem>2G] span[hosts=1] rusage[mem=2G]' -a 'docker(kboltonlab/annotate_wes_ch:3.1)' -J "ann_${folder}[2001-3000]" /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/annotate_pd_scratch_SAMPLE.sh ${folder}
bsub -o logs/final_${folder}_%J_%I.out -G compute-timley -g /bwileytest2 -q general -M 2G -R 'select[mem>2G] span[hosts=1] rusage[mem=2G]' -a 'docker(kboltonlab/annotate_wes_ch:3.1)' -J "ann_${folder}[3001-4000]" /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/annotate_pd_scratch_SAMPLE.sh ${folder}
bsub -o logs/final_${folder}_%J_%I.out -G compute-timley -g /bwileytest2 -q general -M 2G -R 'select[mem>2G] span[hosts=1] rusage[mem=2G]' -a 'docker(kboltonlab/annotate_wes_ch:3.1)' -J "ann_${folder}[4001-4050]" /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/annotate_pd_scratch_SAMPLE.sh ${folder}


for dir in $(ls /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/*/merged/merged.final.tsv | cut -d/ -f12 | grep -v 10); do
    cd /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/RERUN/$dir/merged/
    awk -F'\t' '{print $189,$190,$202,$203,$204,$205,$206,$207}' OFS='\t' merged.final.tsv > merged.annotation.tsv
done
COL merged.final.tsv "key\|gene_aachange\|COSMIC_ID\|CosmicCount\|heme_cosmic_count\|myeloid_cosmic_count\|oncoKB\|sOncogenic"

awk -F'\t' '{print $186,$187,$199,$200,$201,$202,$203,$204}' OFS='\t' merged.final.tsv > merged.annotation.tsv


for dir in {10..11}; do
    echo $dir
    cd $EXOME/results/RERUN/$dir/merged/combined
    ../../../hs.sh
    echo
done