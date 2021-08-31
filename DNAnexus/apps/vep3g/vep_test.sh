/opt/vep/src/ensembl-vep/vep \
    --format vcf \
    -i examples/homo_sapiens_GRCh38.vcf \
    --fork 4 \
    --terms SO \
    --transcript_version \
    --offline \
    --cache \
    --symbol \
    --vcf \
    -o homo_sapiens_GRCh38.vep.vcf \
    --dir /opt/vep/.vep/ \
    --sift p \
    --polyphen p \
    --coding_only \
    --pick \
    --flag_pick \
    --pick_allele \
    --per_gene \
    --pick_allele_gene \
    --flag_pick_allele \
    --flag_pick_allele_gene \
    --plugin Frameshift \
    --plugin Wildtype \
    --everything 1 \
    --assembly GRCh38 \
    --species homo_sapiens \
    --merged \
    --check_existing \
    --force_overwrite && bgzip homo_sapiens_GRCh38.vep.vcf && tabix homo_sapiens_GRCh38.vep.vcf.gz