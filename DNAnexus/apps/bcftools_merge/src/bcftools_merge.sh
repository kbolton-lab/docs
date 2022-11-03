#!/bin/bash
# bcftools_isec_complement 0.0.1

main() {
    set -ex -o pipefail

    echo "Value of reference: '$reference'"
    echo "Value of reference_index_path: '$reference_index_path'"
    echo "Value of tumor_sample_name: '$tumor_sample_name'"
    echo "Value of vcfs: '${vcfs[@]}'"
    echo "Value of dockerimage_bcftools: '$dockerimage_bcftools'"

    # /usr/bin/bcftools -v
    #   "execDepends": [{"name": "bcftools"},
    #   {"name": "tabix"}]
    
    # dx run applet-GFykVk0JQ2896Kzp4bj1707b -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/bcftools_merge/bcftools_merge.json --destination $KELLY:/Brian/test/ -y

   
    dx-download-all-inputs --parallel
    # mv $reference_path .
    # mv $reference_index_path .
    dx download "$reference"
    dx download "$reference_index"


    vcf_string=""
    for vcf in "${vcfs_path[@]}"; do
        vcf_filename=$(basename -- "$vcf")
        mv $vcf .
        # extension="${vcf_filename##*.}"
        # echo $extension
        # if [[ $extension == "gz" ]]; then
        #     vcf_string="$vcf_string $vcf_filename" 
        # fi   
    done

    ls *.vcf.gz
    ls $reference_name
    ls "${reference_name}.fai"
  
    docker load -i $dockerimage_bcftools_path
    #-f $reference_path # 1701353
    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -w /home/dnanexus kboltonlab/sam_bcftools_tabix_bgzip:1.0 \
        /bin/bash -c "/usr/local/bin/bcftools concat --allow-overlaps --remove-duplicates --threads 8 *.vcf.gz | /usr/local/bin/bcftools view -s ${tumor_sample_name} --threads 8 | bcftools filter -i 'REF==\"A\" || REF==\"C\" || REF==\"G\" || REF==\"T\" || REF==\"N\"' | /usr/local/bin/bcftools norm -f $reference_name -m -any --threads 8 -Oz -o $tumor_sample_name.$output_vcf_suffix; tabix $tumor_sample_name.$output_vcf_suffix"
        

    merged_vcf=$(dx upload "$tumor_sample_name.$output_vcf_suffix" --brief)
    dx-jobutil-add-output merged_vcf --class=file "$merged_vcf" 

    merged_vcf_index=$(dx upload "$tumor_sample_name.$output_vcf_suffix.tbi" --brief)
    dx-jobutil-add-output merged_vcf_index --class=file "$merged_vcf_index"
}
