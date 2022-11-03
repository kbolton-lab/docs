#!/bin/bash
# vardict 0.0.1

main() {
    set -ex -o pipefail
    # dx run applet-GFyyvb8JQ2869vky4bJ260xv -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/varscan/vardict_test_inputs_WGS_25.json -y --destination $KELLY:/Brian/test/
    if [[ ! -f /usr/local/bin/Varscan.sh ]]; then
        echo 'File "/usr/local/bin/Varscan.sh" is not there, aborting.'
        exit
    fi
    if [[ ! -f /varscan/VarScan.v2.4.2.jar ]]; then
        echo 'File "/varscan/VarScan.v2.4.2.jar" is not there, aborting.'
        exit
    fi
    which bcftools
    which bgzip
    which tabix

    echo "Value of reference: '${reference}'"
    echo "Value of reference_index: '${reference_index}'"
    echo "Value of reference_dict: '${reference_dict}'"
    echo "Value of array_of_scattered_input: '${interval_beds_array[@]}'" 
    echo "Value of dockerimage_varscan: '${dockerimage_varscan}'" 
    echo "Value of dockerimage_bcftools: '${dockerimage_bcftools}'"

    #dx download "$dockerimage_varscan"
    #ls
    #docker load -i ${dockerimage_varscan_path}
    #docker load -i ${dockerimage_varscan_name}

    # docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin/:/usr/local/bin -v /varscan:/varscan  -w /home/dnanexus kboltonlab/varscan2:1.0
    #/bin/bash /usr/local/bin/Varscan2.sh

    dx download "$tumor_bam"
    dx download "$normal_bam"
    
    tumor_sample_name=$(/usr/bin/samtools view -H $tumor_bam_name | /usr/bin/perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)
    echo $tumor_sample_name
    
    normal_sample_name=$(/usr/bin/samtools view -H $normal_bam_name | /usr/bin/perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1)
    echo $normal_sample_name
   
    interval_files=()
    for file in "${interval_beds_array[@]}"; do
        interval_files+=("-iarray_of_scattered_input=${file}")
    done

 
    map_job=$(dx-jobutil-new-job map \
        "${interval_files[@]}" \
        -i "reference=${reference}" \
        -i "reference_index=${reference_index}" \
        -i "reference_dict=${reference_dict}" \
        -i "tumor_bam=${tumor_bam}" \
        -i "tumor_bai=${tumor_bai}" \
        -i "normal_bam=${normal_bam}" \
        -i "normal_bai=${normal_bai}" \
        -i "dockerimage_varscan=${dockerimage_varscan}" \
        -i "AF_THR=${AF_THR}" \
        -i "tumor_sample_name=${tumor_sample_name}" \
        -i "normal_sample_name=${normal_sample_name}")

    ## add tumor_sample_name for extracting without need to create new instance for "bcftools_extract_tumor" app
    ## get <eid>_23153_0_0 from <eid>_23153_0_0.bqsr.bam
    eid_nameroot=$(echo $tumor_bam_name | cut -d'.' -f1)
    output_name="$eid_nameroot.merged.varcan.vcf.gz"
    postprocess_job=$(dx-jobutil-new-job postprocess \
        -i process_outputs:array:jobref=${map_job}:process_outputs \
        -i tumor_sample_name:string=${tumor_sample_name} \
        -i output_name:string=${output_name} \
        -i "reference:file=${reference}" \
        -i "reference_index:file=${reference_index}" \
        -i "dockerimage_bcftools:file=${dockerimage_bcftools}" \
        --depends-on $map_job)


    dx-jobutil-add-output vcf --class=jobref "$postprocess_job":vcf 
    dx-jobutil-add-output vcf_index --class=jobref "$postprocess_job":vcf_index
    dx-jobutil-add-output tumor_sample_name --class=string "$tumor_sample_name"
}

map() {
    set -ex -o pipefail
    
    echo "Value of array_of_scattered_input: '${array_of_scattered_input[@]}'" 
    echo "Value of dockerimage_varscan: '${dockerimage_varscan}'" 
    echo "Value of AF_THR: '${AF_THR}'" 
    echo "Value of tumor_sample_name: '${tumor_sample_name}'" 
    echo "Value of normal_sample_name: '${normal_sample_name}'" 

    ## might need counter i.e. a "shard" like Cromwell
    process_jobs=()
    i=1
    for interval_file in "${array_of_scattered_input[@]}"
    do
        echo "Value of interval_file: '${interval_file}'"
        
        process_jobs+=($(dx-jobutil-new-job process -i "interval_file=${interval_file}" -i "reference=${reference}" -i "reference_index=${reference_index}" -i "reference_dict=${reference_dict}" -i "tumor_bam=${tumor_bam}" -i "tumor_bai=${tumor_bai}" -i "normal_bam=${normal_bam}" -i "normal_bai=${normal_bai}" -i shard:string=shard-$i -i "dockerimage_varscan=${dockerimage_varscan}" -i "AF_THR=${AF_THR}" -i "tumor_sample_name=${tumor_sample_name}" -i "normal_sample_name=${normal_sample_name}"))
      
        echo "${process_jobs[@]}"
        i=$((i+1))
        echo "$i.after"
    done

    for process_job in "${process_jobs[@]}"
    do
        dx-jobutil-add-output process_outputs --class=array:jobref "${process_job}:vcf" --array
        dx-jobutil-add-output process_outputs --class=array:jobref "${process_job}:vcf_index" --array 
    done

}

process() {
    set -ex -o pipefail

    echo "Value of AF_THR: '${AF_THR}'" 
    echo "Value of tumor_sample_name: '${tumor_sample_name}'" 
    echo "Value of normal_sample_name: '${normal_sample_name}'" 
    echo "Value of interval_file_name: '${interval_file_name}'"
    echo "Value of dockerimage_varscan_path: '${dockerimage_varscan_path}'"

    dx-download-all-inputs --parallel
    mv $tumor_bai_path ~/in/tumor_bam
    mv $normal_bai_path ~/in/normal_bam
    mv $reference_index_path ~/in/reference
    mv $reference_dict_path ~/in/reference

    docker load -i ${dockerimage_varscan_path}

    # docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin/:/usr/local/bin -w /home/dnanexus kboltonlab/vardict:bedtools \
    #docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -v /usr/local/bin/:/usr/local/bin -w /home/dnanexus kboltonlab/varscan2
    /bin/bash /usr/local/bin/Varscan.sh ${reference_path} ${AF_THR} ${tumor_bam_path} ${tumor_sample_name} ${interval_file_path} ${normal_bam_path} ${normal_sample_name} ${shard}.varscan ${interval_file_name}

    vcf=$(dx upload ${shard}.varscan.concat.reheader.vcf.gz --brief)
    dx-jobutil-add-output vcf --class=file "$vcf" 
    dx upload ${shard}.varscan.concat.reheader.vcf.gz --destination project-G4qpk1jJQ285yvbXPFZKXkk8:/Brian/test/varscan/
    
    vcf_index=$(dx upload ${shard}.varscan.concat.reheader.vcf.gz.tbi --brief)
    dx-jobutil-add-output vcf_index --class=file "$vcf_index"
    dx upload ${shard}.varscan.concat.reheader.vcf.gz.tbi --destination project-G4qpk1jJQ285yvbXPFZKXkk8:/Brian/test/varscan/
}

postprocess() {
    set -ex -o pipefail
    if [[ -f /usr/bin/bcftools && -f /usr/bin/tabix ]]; then
        /usr/bin/bcftools -v
        /usr/bin/tabix --version
    else 
        exit
    fi
    echo "Value of process_outputs: '${process_outputs[@]}'"
    echo "Value of output_name: '${output_name}'"
    echo "Value of tumor_sample_name: '${tumor_sample_name}'"
    echo "Value of reference: '${reference}'"
    echo "Value of reference_index: '${reference_index}'"
    echo "Value of dockerimage_bcftools: '${dockerimage_bcftools}'"
    
    dx download "$reference"
    dx download "$reference_index"
    dx download "$dockerimage_bcftools"
    /usr/bin/bcftools -v

    # download vcfs and vcf indices
    mkdir -p /tmp/varscan
    cd /tmp/varscan
    for process_output in "${process_outputs[@]}" 
    do
        dx download "$process_output"
        echo "$process_output_path"
    done
    ls -1sh /tmp/varscan/*.vcf.gz*

    which zcat
    for file in $(ls *.vcf.gz); do
        echo $file
        zcat $file | grep -v "^#" | wc -l
        echo
    done

    cd /home/dnanexus
    
    # extract tumor before bcbio filter or both normal and tumor will be canidates for the filter
    # /usr/bin/bcftools concat --allow-overlaps --remove-duplicates --threads 8 /tmp/vardict/*.vcf.gz | /usr/bin/bcftools view -s ${tumor_sample_name} --threads 8 | bcftools norm -f $reference_name -m -any --threads 8 | /usr/bin/bcftools filter -i "((FMT/AF * FMT/DP < 6) && ((FMT/MQ < 55.0 && FMT/NM > 1.0) || (FMT/MQ < 60.0 && FMT/NM > 2.0) || (FMT/DP < 10) || (FMT/QUAL < 45)))" -m+ -s "BCBIO" --threads 8 -Oz -o $output_name && /usr/bin/tabix $output_name 
    ls $reference_name
    ls "${reference_name}.fai"
  
    docker load -i $dockerimage_bcftools_name
    #-f $reference_path
    # docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -w /home/dnanexus kboltonlab/sam_bcftools_tabix_bgzip:1.0 \
    #     /bin/bash -c "bcftools concat --allow-overlaps --remove-duplicates --threads 8 *.vcf.gz | bcftools view -s ${tumor_sample_name} --threads 8 | bcftools norm -f $reference_name -m -any --threads 8 -Oz -o $output_name && tabix $output_name"
    docker run --rm -v /tmp/varscan:/tmp/varscan -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -w /home/dnanexus kboltonlab/sam_bcftools_tabix_bgzip:1.0 \
        /bin/bash -c "bcftools concat --allow-overlaps --remove-duplicates --threads 8 /tmp/varscan/*.vcf.gz | bcftools view -s ${tumor_sample_name} --threads 8 | bcftools filter -i 'REF==\"A\" || REF==\"C\" || REF==\"G\" || REF==\"T\" || REF==\"N\"' --threads 8 | bcftools norm -f $reference_name -m -any --threads 8 -Oz -o $output_name && tabix $output_name"
    # extract tumor
    # /usr/bin/bcftools concat --allow-overlaps --remove-duplicates --threads 8 /tmp/varscan/*.vcf.gz | /usr/bin/bcftools view -s ${tumor_sample_name} --threads 8 | bcftools norm -f $reference_name -m -any --threads 8 -Oz -o $output_name && /usr/bin/tabix $output_name

    vcf=$(dx upload "$output_name" --brief)
    dx-jobutil-add-output vcf --class=file "$vcf" 

    vcf_index=$(dx upload "$output_name.tbi" --brief)
    dx-jobutil-add-output vcf_index --class=file "$vcf_index"

}
