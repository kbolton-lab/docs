#!/bin/bash
# vardict multithreaded 0.0.1

main() {
    set -ex -o pipefail

    echo "Value of bed: '${bed}'"
    echo "Value of tumor_bam: '${tumor_bam_name}'"
    echo "Value of normal_bam: '${normal_bam_name}'"
    echo "Value of dockerimage_vardict: '${dockerimage_vardict}'"
    echo "Value of AF_THR: '${AF_THR}'" 

    dx-download-all-inputs --parallel

    mv $tumor_bai_path ~/in/tumor_bam
    mv $normal_bai_path ~/in/normal_bam
    mv $reference_index_path ~/in/reference
    mv $reference_dict_path ~/in/reference

    ## add tumor_sample_name for extracting without need to create new instance for "bcftools_extract_tumor" app
    ## TCGA-blah-blah-.bqsr.bam -> TCGA-blah-blah-.vardict.BCBIOfiltered.vcf.gz
    eid_nameroot=$tumor_bam_prefix
    docker load -i ${dockerimage_vardict_path}

    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin/:/usr/local/bin -w /home/dnanexus kboltonlab/vardictjava:1.0 \
        /bin/bash /usr/local/bin/VardictJava_MT.sh ${reference_path} ${AF_THR} ${tumor_bam_path} ${bed_path} ${normal_bam_path} $eid_nameroot


    vcf=$(dx upload "$eid_nameroot.vardict.BCBIOfiltered.vcf.gz" --brief)
    vcf_index=$(dx upload "$eid_nameroot.vardict.BCBIOfiltered.vcf.gz.tbi" --brief)
   
    dx-jobutil-add-output vcf --class=file "$vcf" 
    dx-jobutil-add-output vcf_index --class=file "$vcf_index"
}

