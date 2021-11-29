#!/bin/bash
# alex for archer 0.0.1


main() {
    set -exou pipefail
    wget https://raw.githubusercontent.com/griffithlab/pVACtools/master/tools/pvacseq/VEP_plugins/Wildtype.pm
    echo "Value of vcf: '$vcf'"
    ## the mutect vcf
    echo "Value of vcf: '$intersect_vcf'"
    echo "Value of caller: '$caller'"


    dx-download-all-inputs --parallel
    mv $vcf_index_path ~/in/vcf  
    mv $intersect_vcf_index_path ~/in/intersect_vcf  
    mv $reference_index_path ~/in/reference
    mv $gnomad_file_index_path ~/in/gnomad_file
    mv $gnomad_V3_index_path ~/in/gnomad_V3
    mv $clinvar_file_index_path ~/in/clinvar_file

    ## get <eid>_23153_0_0 from <eid>_23153_0_0.<suffix>
    eid_nameroot=$(echo $vcf_name | cut -d'.' -f1)
    export name="$eid_nameroot"."$caller".final.annotated.vcf.gz

    /usr/bin/bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%AD]\n" $intersect_vcf_path > mutect.tsv
    /usr/bin/bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%AD]\t[%MQ]\n" $vcf_path > vardict.tsv
    
    ## add vardict intersect mutect step
    ## need to add Alex's script at end
    docker load -i ${dockerimage_julia_path}

    # intersect_vcf_path = query or mutect
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/helper:/usr/local/helper -w /home/dnanexus kboltonlab/julia:brian /usr/local/bin/julia /usr/local/helper/SQO.jl mutect.tsv vardict.tsv results.tsv

    sed -i 's/PooledArrays.PooledVector{String, UInt32, Vector{UInt32}}//g' results.tsv
    awk -F'\t' '$13==1 && $7 >= 55 {print}' results.tsv > results.filtered.tsv

    # Annotate
    # Vardict is the query and mutect is the subject
    printf "##INFO=<ID=calpos,Number=1,Type=String,Description=\"Positions of subject vcf for complex variants\">\n##INFO=<ID=calref,Number=1,Type=String,Description=\"REF for Positions of subject vcf for complex variants\">\n##INFO=<ID=calalt,Number=1,Type=String,Description=\"ALT for Positions of subject vcf for complex variants\">\n##INFO=<ID=calpos_best,Number=1,Type=String,Description=\"Best Positions of subject vcf for complex variants\">\n##INFO=<ID=calref_best,Number=1,Type=String,Description=\"Best REF for Positions of subject vcf for complex variants\">\n##INFO=<ID=calalt_best,Number=1,Type=String,Description=\"Best ALT for Positions of subject vcf for complex variants\"" > complex.header;

    /usr/bin/bgzip -f results.filtered.tsv
    /usr/bin/tabix -f -s1 -b2 -e2 results.filtered.tsv.gz
    /usr/bin/bcftools annotate --threads 4 -a results.filtered.tsv.gz -h complex.header -c CHROM,POS,REF,ALT,-,-,-,calpos,calref,calalt,-,-,-,calpos_best,calalt_best,calref_best $vcf_path -Oz -o complex.variant.query.vcf.gz && /usr/bin/tabix complex.variant.query.vcf.gz
    bcftools filter -i 'INFO/calpos~".*"' complex.variant.query.vcf.gz -Oz -o complex.only.variant.query.vcf.gz && /usr/bin/tabix complex.only.variant.query.vcf.gz

    docker load -i $dockerimage_vep_path
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -w /home/dnanexus kboltonlab/vep \
        /bin/bash -c "/opt/vep/src/ensembl-vep/vep \
            --format vcf \
            -i complex.only.variant.query.vcf.gz \
            --fork 4 \
            --terms SO \
            --transcript_version \
            --symbol \
            --vcf \
            --database \
            -o complex.only.variant.query.vep.annotated.vcf \
            --fasta $reference_path \
            --dir /opt/vep/.vep/ \
            --synonyms $synonyms_path \
            --sift p \
            --polyphen p \
            --coding_only \
            --pick \
            --plugin Frameshift \
            --plugin Wildtype \
            --everything 1 \
            --assembly GRCh38 \
            --species homo_sapiens \
            --merged \
            --check_existing \
            --buffer_size 1000 \
            --af_gnomad \
            --custom $gnomad_file_path,gnomADe,vcf,exact,1,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS \
            --custom $gnomad_V3_path,gnomADg,vcf,exact,1,AF,AF_ami,AF_oth,AF_afr,AF_sas,AF_asj,AF_fin,AF_amr,AF_nfe,AF_eas \
            --custom $clinvar_file_path,clinvar,vcf,exact,1,CLINSIGN,PHENOTYPE,SCORE,RCVACC,TESTEDINGTR,PHENOTYPELIST,NUMSUBMIT,GUIDELINES \
            --force_overwrite && bgzip complex.only.variant.query.vep.annotated.vcf && tabix complex.only.variant.query.vep.annotated.vcf.gz"


    /usr/bin/bcftools isec -n+2 -w1 complex.variant.query.vcf.gz ${intersect_vcf_path} -Oz -o $name && /usr/bin/tabix $name

    annotated_vcf=$(dx upload $name --brief)
    annotated_vcf_index=$(dx upload $name.tbi --brief)

    dx-jobutil-add-output annotated_vcf "$annotated_vcf" --class=file
    dx-jobutil-add-output annotated_vcf_index "$annotated_vcf_index" --class=file
}
