#!/bin/bash

mutect=$1
vardict=$2
out=$3

/usr/local/bin/bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%AD]\n" $mutect > mutect.tsv
/usr/local/bin/bcftools query -f "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%AD]\n" $vardict > vardict.tsv

/usr/local/bin/julia /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/vardict_pon2at2percent/resources/usr/local/helper/SQO.jl mutect.tsv vardict.tsv $out

# docker run -it --rm -v /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/test:/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/test -v /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/vardict_pon2at2percent/resources/usr/local/helper:/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/vardict_pon2at2percent/resources/usr/local/helper kboltonlab/bst:1.1 /bin/bash /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/vardict_pon2at2percent/resources/usr/local/helper/helper.sh /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/test/5977335_23153_0_0.mapq0.soft-filtered.gz /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/test/5977335_23153_0_0.vardict.BCBIOfiltered.2.vcf.gz results.tsv