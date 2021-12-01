set -o pipefail
set -o errexit

/gatk/gatk HaplotypeCaller --java-options "-Xmx32g" \
    -R $1 \
    -I $2 \
    -O $3 \
    -L $4 \
    -ERC GVCF \
    -G AS_StandardAnnotation \
    -G StandardAnnotation