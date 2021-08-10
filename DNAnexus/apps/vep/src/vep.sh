#!/bin/bash
# vep 0.0.1
# Generated by dx-app-wizard.
#
# Basic execution pattern: Your app will run on a single machine from
# beginning to end.
#
# Your job's input variables (if any) will be loaded as environment
# variables before this script runs.  Any array inputs will be loaded
# as bash arrays.
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.
#
# See https://documentation.dnanexus.com/developer for tutorials on how
# to modify this file.

main() {

    echo "Value of vcf: '$vcf'"
    echo "Value of vcf_index: '$vcf_index'"
    echo "Value of reference: '$reference'"
    echo "Value of reference_index: '$reference_index'"
    echo "Value of synonyms: '$synonyms'"
    echo "Value of gnomad_file: '$gnomad_file'"
    echo "Value of gnomad_file_index: '$gnomad_file_index'"
    echo "Value of clinvar_file: '$clinvar_file'"
    echo "Value of clinvar_file_index: '$clinvar_file_index'"
    echo "Value of species: '$species'"

    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".

    dx download "$vcf"
    dx download "$vcf_index"
    dx download "$reference"
    dx download "$reference_index"

    if [ -n "$synonyms" ]
    then
        dx download "$synonyms"
    fi
    if [ -n "$gnomad_file" ]
    then
        dx download "$gnomad_file"
    fi
    if [ -n "$gnomad_file_index" ]
    then
        dx download "$gnomad_file_index"
    fi
    if [ -n "$clinvar_file" ]
    then
        dx download "$clinvar_file"
    fi
    if [ -n "$clinvar_file_index" ]
    then
        dx download "$clinvar_file_index"
    fi

    # Fill in your application code here.
    #
    # To report any recognized errors in the correct format in
    # $HOME/job_error.json and exit this script, you can use the
    # dx-jobutil-report-error utility as follows:
    #
    #   dx-jobutil-report-error "My error message"
    #
    # Note however that this entire bash script is executed with -e
    # when running in the cloud, so any line which returns a nonzero
    # exit code will prematurely exit the script; if no error was
    # reported in the job_error.json file, then the failure reason
    # will be AppInternalError with a generic error message.
    #dx load -i /vep.tar.gz

    docker run --rm -v /home/dnanexus:/home/dnanexus -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 -w /home/dnanexus kboltonlab/vep2 \
        /bin/bash -c "/opt/vep/src/ensembl-vep/vep \
            --format vcf \
            -i $vcf_name \
            --fork 4 \
            --terms SO \
            --transcript_version \
            --offline \
            --cache \
            --symbol \
            --vcf \
            -o $vcf_prefix.annotated.vcf \
            --fasta $reference_name \
            --dir /opt/vep/.vep/ \
            --synonyms $synonyms_name \
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
            --custom $gnomad_file_name,gnomADe,vcf,exact,1,AF,AF_AFR,AF_AMR,AF_ASJ,AF_EAS,AF_FIN,AF_NFE,AF_OTH,AF_SAS \
            --custom $clinvar_file_name,clinvar,vcf,exact,1,CLINSIGN,PHENOTYPE,SCORE,RCVACC,TESTEDINGTR,PHENOTYPELIST,NUMSUBMIT,GUIDELINES \
            --force_overwrite && bgzip $vcf_prefix.annotated.vcf && tabix $vcf_prefix.annotated.vcf.gz"

    # The following line(s) use the dx command-line tool to upload your file
    # outputs after you have created them on the local file system.  It assumes
    # that you have used the output field name for the filename for each output,
    # but you can change that behavior to suit your needs.  Run "dx upload -h"
    # to see more options to set metadata.

    annotated_vcf=$(dx upload  $vcf_prefix.annotated.vcf.gz --brief)
    annotated_vcf_index=$(dx upload  $vcf_prefix.annotated.vcf.gz.tbi --brief)

    # The following line(s) use the utility dx-jobutil-add-output to format and
    # add output variables to your job's output as appropriate for the output
    # class.  Run "dx-jobutil-add-output -h" for more information on what it
    # does.

    dx-jobutil-add-output annotated_vcf "$annotated_vcf" --class=file
    dx-jobutil-add-output annotated_vcf_index "$annotated_vcf_index" --class=file
}