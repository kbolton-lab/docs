#!/bin/bash
# bcftools_isec_complement 0.0.1

main() {
    set -ex -o pipefail
    
    echo "Value of vcf_in: '$vcf_in'"
    echo "Value of vcf_in_index: '$vcf_in_index'"
    echo "Value of exclude_vcf: '$exclude_vcf'"
    echo "Value of exclude_vcf_index: '$exclude_vcf_index'"
    echo "Value of output_type: '$output_type'"
    echo "Value of output_vcf_name: '$output_vcf_name'"

   
    dx-download-all-inputs --parallel
    echo $vcf_in_index_path
    mv $vcf_in_index_path ~/in/vcf_in
    echo "ls ~/in/vcf_in"
    ls ~/in/vcf_in
    echo
    
    mv $exclude_vcf_index_path ~/in/exclude_vcf
    echo "ls ~/in/exclude_vcf"
    ls ~/in/exclude_vcf
  
    docker load -i /bcftools.tar.gz
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -w /home/dnanexus kboltonlab/sam_bcftools_tabix_bgzip:1.0 bash -c "/usr/local/bin/bcftools isec -C -w1 ${vcf_in_path} ${exclude_vcf_path} -O${output_type} -o ${output_vcf_name} && tabix ${output_vcf_name}" 

    complement_vcf=$(dx upload ${output_vcf_name} --brief)
    complement_vcf_index=$(dx upload ${output_vcf_name}.tbi --brief)

    dx-jobutil-add-output complement_vcf "$complement_vcf" --class=file
    dx-jobutil-add-output complement_vcf_index "$complement_vcf_index" --class=file
}
