#!/bin/bash
# bcftools_isec_complement 0.0.1

main() {
    set -ex -o pipefail
    
    echo "Value of vcf1: '$vcf1'"
    echo "Value of vcf1_index: '$vcf1_index'"
    echo "Value of vcf2: '$vcf2'"
    echo "Value of vcf2_index: '$vcf2_index'"
    echo "Value of output_type: '$output_type'"
    echo "Value of output_vcf_name: '$output_vcf_name'"

   
    dx-download-all-inputs --parallel
    mv $vcf1_index_path ~/in/vcf1
    mv $vcf2_index_path ~/in/vcf2
  
    /usr/bin/bcftools isec -n+2 -w1 ${vcf1_path} ${vcf2_path} -O${output_type} -o ${output_vcf_name} --threads 8 && tabix ${output_vcf_name}

    intersecting_vcf=$(dx upload ${output_vcf_name} --brief)
    intersecting_vcf_index=$(dx upload ${output_vcf_name}.tbi --brief)

    dx-jobutil-add-output intersecting_vcf "$intersecting_vcf" --class=file
    dx-jobutil-add-output intersecting_vcf_index "$intersecting_vcf_index" --class=file
}
