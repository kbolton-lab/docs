#!/bin/bash
# pon2at2percent 0.0.1


main() {
    set -eou pipefail
    echo "Value of vcf: '$vcf'"
    #echo "Value of vcf_index: '$vcf_index'"
    echo "Value of vcf2PON: '$vcf2PON'"
    echo "Value of vcf2PON_index: '$vcf2PON_index'"
    echo "Value of caller: '$caller'"
    echo "Value of tumor_sample_name: '$tumor_sample_name'"


    dx-download-all-inputs --parallel
    #dx-download-all-inputs
    mv $vcf2PON_index_path ~/in/vcf2PON   
  
    export name="$caller"."$tumor_sample_name".final.annotated.vcf.gz

    printf "##INFO=<ID=PON_2AT2_percent,Number=1,Type=Integer,Description=\"If 2 PoN samples have variant at >=2 percent\">\n" > pon2.header;
    /usr/bin/bcftools query -f "%CHROM\t%POS\t%REF\t%REF\t1\n" $vcf2PON_path > normal2.txt
    /usr/bin/bgzip -f normal2.txt
    /usr/bin/tabix -f -s1 -b2 -e2 normal2.txt.gz
    /usr/bin/bcftools annotate -a normal2.txt.gz -h pon2.header -c CHROM,POS,REF,ALT,PON_2AT2_percent $vcf_path -Oz -o $name
    tabix $name


    annotated_pon2_vcf=$(dx upload $name --brief)
    annotated_pon2_vcf_index=$(dx upload $name.tbi --brief)

    dx-jobutil-add-output annotated_pon2_vcf "$annotated_pon2_vcf" --class=file
    dx-jobutil-add-output annotated_pon2_vcf_index "$annotated_pon2_vcf_index" --class=file
}
