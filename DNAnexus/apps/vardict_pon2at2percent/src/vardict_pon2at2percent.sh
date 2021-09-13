#!/bin/bash
# pon2at2percent 0.0.1


main() {
    set -eou pipefail
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
    
    ## add vardict intersect mutect step
    ## need to add Alex's script at end
    /usr/bin/bcftools isec -n+2 -w1 ${vcf_path} ${intersect_vcf_path} -Oz -o $eid_nameroot.vardict.intersect.mutect.vcf.gz
    

    printf "##INFO=<ID=PON_2AT2_percent,Number=1,Type=Integer,Description=\"If 2 PoN samples have variant at >=2 percent\">\n" > pon2.header;
    /usr/bin/bcftools query -f "%CHROM\t%POS\t%REF\t%REF\t1\n" $vcf2PON_path > normal2.txt
    /usr/bin/bgzip -f normal2.txt
    /usr/bin/tabix -f -s1 -b2 -e2 normal2.txt.gz
    /usr/bin/bcftools annotate --threads 4 -a normal2.txt.gz -h pon2.header -c CHROM,POS,REF,ALT,PON_2AT2_percent $eid_nameroot.vardict.intersect.mutect.vcf.gz -Oz -o $name
    tabix $name


    annotated_pon2_vcf=$(dx upload $name --brief)
    annotated_pon2_vcf_index=$(dx upload $name.tbi --brief)

    dx-jobutil-add-output annotated_pon2_vcf "$annotated_pon2_vcf" --class=file
    dx-jobutil-add-output annotated_pon2_vcf_index "$annotated_pon2_vcf_index" --class=file
}
