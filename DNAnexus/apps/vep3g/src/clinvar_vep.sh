#!/bin/bash

 kboltonlab/vep3
clinvar_file_path=/storage1/fs1/bolton/Active/data/hg38/vcf/clinvar_20220816.vcf.gz
clinvar_file_path=/storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/vcf/0001.vcf.gz
/opt/vep/src/ensembl-vep/vep \
            --format vcf \
            -i $clinvar_file_path \
            --fork 4 \
            --terms SO \
            --transcript_version \
            --offline \
            --cache \
            --symbol \
            --vcf \
            --dir_plugins /opt/vep/.vep/Plugins/ \
            -o /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/vcf/0001.annotated.vcf  \
            --fasta $EXOME/GRCh38_full_analysis_set_plus_decoy_hla.fa \
            --dir /storage1/fs1/bolton/Active/pvcfs/jie_test/vep/hg38/cache/ \
            --synonyms $GMSROOT/gc2560/core/model_data/2887491634/build50f99e75d14340ffb5b7d21b03887637/chromAlias.ensembl.txt \
            --plugin Frameshift \
            --plugin Wildtype \
            --assembly GRCh38 \
            --species homo_sapiens \
            --sift p --polyphen p --pick --pick_order canonical,rank,mane,ccds,appris,tsl,biotype,length --everything 1 --merged --check_existing --buffer_size 1000 --af_gnomad \
            --custom $clinvar_file_path,clinvar,vcf,exact,1,CLNSIG,ORIGIN \
            --force_overwrite
            
            bgzip clinvar_20220816.vep.annotated.vcf 
            tabix clinvar_20220816.vep.annotated.vcf.gz


--dir_cache $EXOME_CALLS/pvcf/vep \
--dir_plugins /opt/vep/.vep/Plugins/ \