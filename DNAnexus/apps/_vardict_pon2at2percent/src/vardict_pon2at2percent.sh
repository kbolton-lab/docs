#!/bin/bash
# pon2at2percent 0.0.1


main() {
    set -exou pipefail
    echo "Value of vcf: '$vcf'"
    ## the mutect vcf
    echo "Value of vcf: '$intersect_vcf'"
    #echo "Value of vcf_index: '$vcf_index'"
    echo "Value of vcf2PON: '$vcf2PON'"
    echo "Value of vcf2PON_index: '$vcf2PON_index'"
    echo "Value of caller: '$caller'"
    # echo "Value of tumor_sample_name: '$tumor_sample_name'"


    dx-download-all-inputs --parallel
    mv $vcf_index_path ~/in/vcf  
    mv $intersect_vcf_index_path ~/in/intersect_vcf  
    mv $vcf2PON_index_path ~/in/vcf2PON

    ## get <eid>_23153_0_0 from <eid>_23153_0_0.<suffix>
    eid_nameroot=$(echo $vcf_name | cut -d'.' -f1)
    export name="$eid_nameroot"."$caller".final.annotated.vcf.gz

    # intersect_vcf_path=/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/test/5977335_23153_0_0.mapq0.soft-filtered.gz
    # vcf_path=/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/test/5977335_23153_0_0.vardict.BCBIOfiltered.2.vcf.gz
    /usr/bin/bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%AD]\n" $intersect_vcf_path > mutect.tsv
    /usr/bin/bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%AD]\t[%MQ]\n" $vcf_path > vardict.tsv
    
    ## add vardict intersect mutect step
    ## need to add Alex's script at end
    docker load -i ${dockerimage_julia_path}

    # intersect_vcf_path = query or mutect
    # docker run -it --rm -v /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/vardict_pon2at2percent/resources/usr/local/test:/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/vardict_pon2at2percent/resources/usr/local/test -v /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/vardict_pon2at2percent/resources/usr/local/helper:/usr/local/helper -w /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/vardict_pon2at2percent/resources/usr/local/test kboltonlab/julia:brian /usr/local/bin/julia /usr/local/helper/SQO.jl mutect.tsv vardict.tsv results.tsv 
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

    /usr/bin/bcftools isec -n+2 -w1 complex.variant.query.vcf.gz ${intersect_vcf_path} -Oz -o $eid_nameroot.vardict.intersect.mutect.vcf.gz && /usr/bin/tabix $eid_nameroot.vardict.intersect.mutect.vcf.gz
    
    printf "##INFO=<ID=PON_2AT2_percent,Number=1,Type=Integer,Description=\"If 2 PoN samples have variant at >=2 percent\">\n" > pon2.header;
    /usr/bin/bcftools query -f "%CHROM\t%POS\t%REF\t%REF\t1\n" $vcf2PON_path > normal2.txt
    /usr/bin/bgzip -f normal2.txt
    /usr/bin/tabix -f -s1 -b2 -e2 normal2.txt.gz
    /usr/bin/bcftools concat -a -D $eid_nameroot.vardict.intersect.mutect.vcf.gz complex.only.variant.query.vcf.gz | /usr/bin/bcftools annotate --threads 4 -a normal2.txt.gz -h pon2.header -c CHROM,POS,REF,ALT,PON_2AT2_percent -Oz -o $name
    # /usr/bin/bcftools annotate --threads 4 -a normal2.txt.gz -h pon2.header -c CHROM,POS,REF,ALT,PON_2AT2_percent $eid_nameroot.vardict.intersect.mutect.vcf.gz -Oz -o $name
    /usr/bin/tabix $name


    annotated_pon2_vcf=$(dx upload $name --brief)
    annotated_pon2_vcf_index=$(dx upload $name.tbi --brief)

    dx-jobutil-add-output annotated_pon2_vcf "$annotated_pon2_vcf" --class=file
    dx-jobutil-add-output annotated_pon2_vcf_index "$annotated_pon2_vcf_index" --class=file
}
