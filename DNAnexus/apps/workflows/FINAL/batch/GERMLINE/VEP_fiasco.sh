bsub -oo test.bsub.log -g /bwileytest3 -G compute-timley -q general -M 16G -R 'select[mem>16G] span[hosts=1] rusage[mem=16G]' -a "docker(kboltonlab/jie_vep:3.0)" /bin/bash -c "/opt/vep/src/ensembl-vep/vep \
   --format vcf \
   -i up_or_down_stream_gene_variant.vcf.gz \
   -o up_or_down_stream_gene_variant.annotated_no_order.vcf \
   --fork 8 --terms SO --transcript_version \
   --offline \
   --cache \
   --symbol \
   --vcf \
   --fasta /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/GRCh38_full_analysis_set_plus_decoy_hla.fa \
   --dir_plugins /opt/vep/.vep/Plugins/ \
   --dir_cache /storage1/fs1/bolton/Active/projects/mocha/UKBB/ukbb_calls/pvcf/vep \
   --synonyms /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build50f99e75d14340ffb5b7d21b03887637/chromAlias.ensembl.txt \
   --plugin Frameshift \
   --plugin Wildtype \
   --assembly GRCh38 \
   --species homo_sapiens \
   --sift p \
   --polyphen p \
   --pick \
   --everything 1 \
   --merged \
   --check_existing \
   --buffer_size 1000 \
   --af_gnomad \
   --custom /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/genome-db-ensembl-gnomad/2dd4b53431674786b760adad60a29273/fixed_b38_exome.vcf.gz,gnomADe,vcf,exact,1,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS \
   --custom /storage1/fs1/bolton/Active/data/hg38/vcf/gnomad.genomes.r3.0.sites.full_trimmed_info.af_only.UKBB_intervals.vcf.bgz,gnomADg,vcf,exact,1,AF,AF_ami,AF_oth,AF_afr,AF_sas,AF_asj,AF_fin,AF_amr,AF_nfe,AF_eas \
   --custom /storage1/fs1/bga/Active/gmsroot/gc2560/core/custom_clinvar_vcf/v20181028/custom.vcf.gz,clinvar,vcf,exact,1,CLINSIGN,PHENOTYPE,SCORE,RCVACC,TESTEDINGTR,PHENOTYPELIST,NUMSUBMIT,GUIDELINES \
   --force_overwrite"

bsub -Is -g /bwileytest3 -G compute-timley -q general-interactive -M 16G -R 'select[mem>16G] span[hosts=1] rusage[mem=16G]' -a "docker(kboltonlab/vep3)" /bin/bash
bsub -Is -g /bwileytest3 -G compute-timley -q general-interactive -M 16G -R 'select[mem>16G] span[hosts=1] rusage[mem=16G]' -a "docker(ensemblorg/ensembl-vep:release_104.3)" /bin/bash
bsub -Is -g /bwileytest3 -G compute-timley -q general-interactive -M 16G -R 'select[mem>16G] span[hosts=1] rusage[mem=16G]' -a "docker(kboltonlab/jie_vep:3.0)" /bin/bash


canonical,mane,rank,ccds,appris,tsl,biotype,length

/opt/vep/src/ensembl-vep/vep \
   --format vcf \
   -i /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/germline/NSAMPS/combined.pon.germline2.stream_gene_variant.vcf \
   -o /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/germline/NSAMPS/combined.pon.germline2.stream_gene_variant.annotated.vcf \
   --fork 8 --terms SO --transcript_version \
   --offline \
   --cache \
   --symbol \
   --vcf \
   --fasta /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/GRCh38_full_analysis_set_plus_decoy_hla.fa \
   --dir_plugins /opt/vep/.vep/Plugins/ \
   --dir_cache /storage1/fs1/bolton/Active/projects/mocha/UKBB/ukbb_calls/pvcf/vep \
   --synonyms /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build50f99e75d14340ffb5b7d21b03887637/chromAlias.ensembl.txt \
   --plugin Frameshift \
   --plugin Wildtype \
   --assembly GRCh38 \
   --species homo_sapiens \
   --sift p \
   --polyphen p \
   --pick \
   --pick_order canonical,mane,rank,ccds,appris,tsl,biotype,length \
   --everything 1 \
   --merged \
   --check_existing \
   --buffer_size 1000 \
   --af_gnomad \
   --custom /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/genome-db-ensembl-gnomad/2dd4b53431674786b760adad60a29273/fixed_b38_exome.vcf.gz,gnomADe,vcf,exact,1,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS \
   --custom /storage1/fs1/bolton/Active/data/hg38/vcf/gnomad.genomes.r3.0.sites.full_trimmed_info.af_only.UKBB_intervals.vcf.bgz,gnomADg,vcf,exact,1,AF,AF_ami,AF_oth,AF_afr,AF_sas,AF_asj,AF_fin,AF_amr,AF_nfe,AF_eas \
   --custom /storage1/fs1/bga/Active/gmsroot/gc2560/core/custom_clinvar_vcf/v20181028/custom.vcf.gz,clinvar,vcf,exact,1,CLINSIGN,PHENOTYPE,SCORE,RCVACC,TESTEDINGTR,PHENOTYPELIST,NUMSUBMIT,GUIDELINES \
   --force_overwrite


/opt/vep/src/ensembl-vep/vep \
   --format vcf \
   -i /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/lung/update/pass/final/up_or_down_stream_gene_variant.vcf.gz \
   -o /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/lung/update/pass/final/up_or_down_stream_gene_variant.annotated_with_order.vcf \
   --fork 8 --terms SO --transcript_version \
   --offline \
   --cache \
   --symbol \
   --vcf \
   --fasta /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/GRCh38_full_analysis_set_plus_decoy_hla.fa \
   --dir_plugins /opt/vep/.vep/Plugins/ \
   --dir_cache /storage1/fs1/bolton/Active/projects/mocha/UKBB/ukbb_calls/pvcf/vep \
   --synonyms /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build50f99e75d14340ffb5b7d21b03887637/chromAlias.ensembl.txt \
   --plugin Frameshift \
   --plugin Wildtype \
   --assembly GRCh38 \
   --species homo_sapiens \
   --sift p \
   --polyphen p \
   --pick \
   --pick_order canonical,mane,rank,ccds,appris,tsl,biotype,length \
   --everything 1 \
   --merged \
   --check_existing \
   --buffer_size 1000 \
   --af_gnomad \
   --custom /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/genome-db-ensembl-gnomad/2dd4b53431674786b760adad60a29273/fixed_b38_exome.vcf.gz,gnomADe,vcf,exact,1,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS \
   --custom /storage1/fs1/bolton/Active/data/hg38/vcf/gnomad.genomes.r3.0.sites.full_trimmed_info.af_only.UKBB_intervals.vcf.bgz,gnomADg,vcf,exact,1,AF,AF_ami,AF_oth,AF_afr,AF_sas,AF_asj,AF_fin,AF_amr,AF_nfe,AF_eas \
   --custom /storage1/fs1/bga/Active/gmsroot/gc2560/core/custom_clinvar_vcf/v20181028/custom.vcf.gz,clinvar,vcf,exact,1,CLINSIGN,PHENOTYPE,SCORE,RCVACC,TESTEDINGTR,PHENOTYPELIST,NUMSUBMIT,GUIDELINES \
   --force_overwrite


original:
--synonyms chromAlias.ensembl.txt \
--coding_only true \
--dir vep_zip
--pick_order canonical,  rank,mane,ccds,  appris,tsl,biotype,length
                         rank=2

_unannotated:
--synonyms chromAlias.ensembl.txt \
 \
--dir vep_zip
--pick_order canonical,  mane,ccds,rank,  appris,tsl,biotype,length 
                                   rank=4 ? wrong?


