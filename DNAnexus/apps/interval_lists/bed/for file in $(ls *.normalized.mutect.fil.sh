for file in $(ls *.normalized.mutect.filtered.vcf.gz | head -n3 | tail -n1); do
   prefix=$(basename $file .normalized.mutect.filtered.vcf.gz)
   bsub -oo logs/$prefix.bcftoolsR.log -g /bwileytest2 -G compute-timley -q general -M 4G -R 'select[mem>4G] span[hosts=1] rusage[mem=4G]' -a "docker(kboltonlab/r_bcftools_vep)" /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/bcftools.sh $prefix $file
done
   tabix -f ../$prefix.mutect.vep.annotated.vcf.gz
   bcftools filter -i 'FILTER="PASS"' $file -Oz -o $prefix.pass.mutect.gz && tabix -f $prefix.pass.mutect.gz && bcftools isec -C -w1 $prefix.pass.mutect.gz $HG38/vcf/af-only-gnomad.biallelic.above.005.leftalign.hg38.UKBB.vcf.gz ../$prefix.mutect.vep.annotated.vcf.gz -Oz -o $prefix.possible.intersect.with.vardict.vcf.gz && tabix $prefix.possible.intersect.with.vardict.vcf.gz
   tabix -f ../$prefix.vardict.final.annotated.vcf.gz
   bcftools isec -w1 -n=2 $prefix.possible.intersect.with.vardict.vcf.gz ../$prefix.vardict.final.annotated.vcf.gz -Oz -o $prefix.intersect.with.vardict.vcf.gz && tabix $prefix.intersect.with.vardict.vcf.gz

   tabix -f ../msk_pileup/$prefix.normal.pileup.vcf.gz
   bcftools +fill-tags -Oz -o $prefix.RD.vcf.gz ../msk_pileup/$prefix.normal.pileup.vcf.gz -- -t "PON_RefDepth=sum(RD)"
   bcftools +fill-tags -Oz -o $prefix.RD_AD.vcf.gz $prefix.RD.vcf.gz -- -t "PON_AltDepth=sum(AD)" && tabix $prefix.RD_AD.vcf.gz
   bcftools annotate --threads 32 -a $prefix.RD_AD.vcf.gz -c PON_RefDepth,PON_AltDepth $prefix.intersect.with.vardict.vcf.gz | bcftools filter -i 'PON_RefDepth>0' | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/PON_RefDepth\t%INFO/PON_AltDepth\t[%AD]\n' > $prefix.fisher.input;
   /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/fisher.R $prefix.fisher.input $prefix.fisher.output
   bgzip -f $prefix.fisher.output
   tabix -f -s1 -b2 -e2 $prefix.fisher.output.gz
   bcftools annotate -a $prefix.fisher.output.gz -h fisher.header -c CHROM,POS,REF,ALT,-,-,-,-,PON_FISHER $prefix.intersect.with.vardict.vcf.gz -Oz -o $prefix.fisherPON.vcf.gz
   printf "TOTAL_VARIANTS: "; bcftools view -H $prefix.fisherPON.vcf.gz | wc -l
   printf "TOTAL_FAILED_ORIG_PON: "; bcftools filter -i "PON_FISHER > 1.260958e-09" $prefix.fisherPON.vcf.gz | bcftools view -H | wc -l

   # clean
   rm $prefix.possible.intersect.with.vardict.vcf.gz* $prefix.intersect.with.vardict.vcf.gz* $prefix.RD.vcf.gz* $prefix.RD_AD.vcf.gz* $prefix.fisher.input;
done


for file in $(ls *.fisherPON.vcf.gz); do
   prefix=$(basename $file.fisherPON.vcf.gz)
   bsub -oo logs/$prefix.vep.log -g /bwileytest2 -G compute-timley -q general -M 8G -R 'select[mem>8G] span[hosts=1] rusage[mem=8G]' -a "docker(kboltonlab/ic_vep:latest)" /bin/bash -c "/usr/bin/perl -I /opt/lib/perl/VEP/Plugins /usr/bin/variant_effect_predictor.pl \
   --format vcf \
   --vcf \
   --fork 4 \
   --terms SO \
   --transcript_version \
   --offline \
   --cache \
   --symbol \
   -o $caller.samples.$sample.gnomAD_AF_filter.filtered.pileup.fisherPON.filtered.pileup.fisherPON_annotated2.vcf \
   -i $caller.samples.$sample.gnomAD_AF_filter.filtered.pileup.fisherPON.filtered.pileup.fisherPON_annotated.vcf.gz \
   --synonyms /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/2887491634/build50f99e75d14340ffb5b7d21b03887637/chromAlias.ensembl.txt \
   --dir /storage1/fs1/bolton/Active/projects/mocha/UKBB/ukbb_calls/pvcf/vep_zip \
   --fasta /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/GCA_000001405.15_GRCh38_full_analysis_set.fna \
   --plugin 'Frameshift' \
   --plugin 'Wildtype' \
   --assembly GRCh38 \
   --cache_version 104 \
   --species homo_sapiens \
   --sift p --polyphen p --pick --pick_order canonical,rank,mane,ccds,appris,tsl,biotype,length --everything 1 --merged --check_existing --buffer_size 1000 --af_gnomad \
   --custom /storage1/fs1/bga/Active/gmsroot/gc2560/core/model_data/genome-db-ensembl-gnomad/2dd4b53431674786b760adad60a29273/fixed_b38_exome.vcf.gz,gnomADe,vcf,exact,1,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS \
   --custom /storage1/fs1/bolton/Active/data/hg38/vcf/gnomad.genomes.r3.0.sites.full_trimmed_info.af_only.above.0.005.vcf.gz,gnomADg,vcf,exact,1,AF,AF_ami,AF_oth,AF_afr,AF_sas,AF_asj,AF_fin,AF_amr,AF_nfe,AF_eas \
   --custom /storage1/fs1/bga/Active/gmsroot/gc2560/core/custom_clinvar_vcf/v20181028/custom.vcf.gz,clinvar,vcf,exact,1,CLINSIGN,PHENOTYPE,SCORE,RCVACC,TESTEDINGTR,PHENOTYPELIST,NUMSUBMIT,GUIDELINES && bgzip -f $caller.samples.$sample.gnomAD_AF_filter.filtered.pileup.fisherPON.filtered.pileup.fisherPON_annotated2.vcf && tabix -f $caller.samples.$sample.gnomAD_AF_filter.filtered.pileup.fisherPON.filtered.pileup.fisherPON_annotated2.vcf.gz"
done

--plugin SpliceAI,snv=/storage1/fs1/bolton/Active/data/hg38/vcf/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/storage1/fs1/bolton/Active/data/hg38/vcf/spliceai_scores.raw.indel.hg38.vcf.gz \







# for file in $(ls *.bqsr.bam); do
#    sample=$(SAMPLE $file)
#    printf "$sample\t$file\n"
# done

# /opt/GetBaseCountsMultiSample/GetBaseCountsMultiSample --fasta /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/GRCh38_full_analysis_set_plus_decoy_hla.fa --bam_fof input.bams.txt --vcf $EXOME/hotspots.bb.cosmic.vcf --output MSK.hotspots.bb.cosmic.vcf --thread 32;
# bgzip -f $msk_out && tabix -f $msk_out.gz

# bcftools +fill-tags -Oz -o RD.vcf.gz $msk_out.gz -- -t "PON_RefDepth=sum(RD)"
# bcftools +fill-tags -Oz -o RD_AD.vcf.gz RD.vcf.gz -- -t "PON_AltDepth=sum(AD)" && tabix RD_AD.vcf.gz