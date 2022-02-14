## Modify anything that is not a function...


```shell
annotate_CH_pd_docker3_Mutect_Vardict_ponChange.UKBB.scratch1.R \
   --variant-calls-files ${file_id}_23153_0_0.mutect.vep.annotated.vcf.gz,${file_id}_23153_0_0.vardict.final.annotated.vcf.gz 
   --truncating /scratch1/fs1/bolton/annotation_files/BB.truncating.more.than.1.tsv \
   --TSG-file /scratch1/fs1/bolton/annotation_files/gene_census_TSG.txt \
   --sample_id $file_id \
   --out_folder /storage1/fs1/bolton/Active/projects/mocha/UKBB/exome/results/$folder/final_outputs \
   --oncoKB-curated /scratch1/fs1/bolton/annotation_files/all_curated_genes_v3.0.tsv \
   --pd-annotation-file /scratch1/fs1/bolton/annotation_files/pd_table_kbreview_bick_trunc3.txt \
   --blacklist /scratch1/fs1/bolton/annotation_files/ENCFF356LFX.bed \
   --segemental-duplications /scratch1/fs1/bolton/annotation_files/dup.grch38.bed.gz \
   --simple-repeats /scratch1/fs1/bolton/annotation_files/simpleRepeat.bed \
   --repeat-masker /scratch1/fs1/bolton/annotation_files/repeatMaskerJoinedCurrent.bed \
   --bolton-bick-vars /scratch1/fs1/bolton/annotation_files/bick.bolton.vars3.txt \
   --mut2-bick /scratch1/fs1/bolton/annotation_files/topmed.n2.mutation.c.p.txt \
   --mut2-kelly /scratch1/fs1/bolton/annotation_files/kelly.n2.mutation.c.p.txt \
   --matches2 /scratch1/fs1/bolton/annotation_files/matches.2.c.p.txt"
   ```
   Reg
   
   ```python
   import python
   ```
