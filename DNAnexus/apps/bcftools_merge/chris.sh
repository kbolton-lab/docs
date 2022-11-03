Docker
1 line 

https://documentation.dnanexus.com/getting-started/tutorials/cli-quickstart

dx-download-all-inputs --parallel

Running of docker:
 #docker load -i $dockerimage_msk_getbasecounts_name
    eid_nameroot=trio1.vcf
    docker load -i $dockerimage_chris_name
    docker run \
        --rm  \
        -v /home/dnanexus:/home/dnanexus \
        -v /mnt/UKBB_Exome_2021:/mnt/UKBB_Exome_2021 \
        -v /usr/local/msk/bin/:/usr/local/msk/bin 
        -w /home/dnanexus \
             kboltonlab/msk_getbasecounts:3.0 \
        /opt/tools/runchris.py -m $reference_name "$bams" $vcf_filename $msk_out -o $eid_nameroot -d path/to/database
    /usr/bin/tabix $eid_nameroot 

    applet script.sh output_vcf_suffix
    # dx upload -h
    vcf_out=$(dx upload $eid_nameroot --brief)
    vcf_out_index=$(dx upload $eid_nameroot.tbi --brief)
    vcf_normal_out=$(dx upload $msk_out.gz --brief)

    # dx-jobutil-add-output -h
    dx-jobutil-add-output vcf_out "$vcf_out" --class=file
    dx-jobutil-add-output vcf_out_index "$vcf_out_index" --class=file
    dx-jobutil-add-output vcf_normal_out "$vcf_normal_out" --class=file



####################################################################################################
IMPORTANT!
dxapp.json

"inputSpec": [
    {
      "name": "bam_mom",
      "label": "bam_mom",
      "class": "file",
      "optional": true,
      "patterns": [
        "*.vcf",
        ".vcf.gz"
      ],
      "help": ""
    },
    {
      "name": "bam_mom_index",
      "label": "bam_mom_index",
      "class": "file",
      "optional": true,
      "patterns": [
        ".vcf.gz.tbi"
      ],
      "help": ""
    },
    {
        "name": "bam_dad",
        "label": "bam_mom",
        "class": "file",
        "optional": true,
        "patterns": [
          "*.vcf",
          ".vcf.gz"
        ],
        "help": ""
      },
      {
        "name": "bam_dad_index",
        "label": "bam_mom_index",
        "class": "file",
        "optional": true,
        "patterns": [
          ".vcf.gz.tbi"
        ],
        "help": ""
      },
      {
        "name": "bam_c",
        "label": "bam_mom",
        "class": "file",
        "optional": true,
        "patterns": [
          "*.vcf",
          ".vcf.gz"
        ],
        "help": ""
      },
      {
        "name": "bam_c_index",
        "label": "bam_mom_index",
        "class": "file",
        "optional": true,
        "patterns": [
          ".vcf.gz.tbi"
        ],
        "help": ""
      },
      {
        "name": "dockerimage_chris",
        "label": "dockerimage_chris",
        "class": "file",
        "optional": false,
        "default": {
            "$dnanexus_link": {
                "project": "project-G3Yj1vjJ6XG579jbKyjXPGGY",
                "id": "file-G4K54K8J6XG9BKKJ6yZyXg9V"
            }
          },
        "patterns": [
          ".vcf.gz.tbi"
        ],
        "help": ""
      },
      {
      "name": "dockerimage_msk_getbasecountsR",
      "class": "file",
      "optional": true,
      "default": {
        "$dnanexus_link": {
            "project": "project-G3Yj1vjJ6XG579jbKyjXPGGY",
            "id": "file-G4K54K8J6XG9BKKJ6yZyXg9V"
        }
      },
      {
      "name": "output_name",
      "label": "output_name",
      "class": "string",
      "optional": true,
      "default": "mutect.filtered.vcf.gz",
      "help": ""
    }
],
"outputSpec": [
    {
      "name": "vcf_normal_out",
      "label": "vcf_normal_out",
      "class": "file",
      "patterns": [
        "*"
      ],
      "help": ""
    },
    {
      "name": "vcf_out",
      "label": "vcf_out",
      "class": "file",
      "patterns": [
        "*.vcf",
        "*.vcf.gz"
      ],
      "help": ""
    },
    {
      "name": "vcf_out_index",
      "label": "vcf_out_index",
      "class": "file",
      "patterns": [
        "*.vcf.gz.tbi"
      ],
      "help": ""
    }
  ]


Machine Specs:
"runSpec": {
    "timeoutPolicy": {
      "*": {
        "hours": 1
      }
    },
    "interpreter": "bash",
    "file": "src/bcftools_isec.sh",
    "execDepends":[
      {"name": "bcftools"},
      {"name": "tabix"}
    ],
    "distribution": "Ubuntu",
    "release": "20.04",
    "version": "0"
  },
  "access": {
    "project": "CONTRIBUTE"
  },
  "regionalOptions": {
    "aws:eu-west-2": {
      "systemRequirements": {
        "*": {
          "instanceType": "mem3_ssd1_v2_x2"
        }
      }
    }
  }


  #############
  ## build app
  dx upload --destination project-G3Yj1vjJ6XG579jbKyjXPGGY:/test/tools/ my.docker.tar.gz
  dx build --destination project-G3Yj1vjJ6XG579jbKyjXPGGY:/test/tools/ annotate_ch_pd -f

  ## run app
  # 1. path to tool
  # 2. input.json
  # 3. --destination folder
  # you can use dx-toolkit CLI 
    # cram=$(dx find data --class file --path UKBB_Exome_2021:"/Bulk/Exome sequences/Exome OQFE CRAM files/$subfolder" --json --name ${sample}_23153_0_0.cram | jq -r '.[].describe.id')
    # cram_index=$(dx find data --class file --path UKBB_Exome_2021:"/Bulk/Exome sequences/Exome OQFE CRAM files/$subfolder" --json --name ${sample}_23153_0_0.cram.crai | jq -r '.[].describe.id')

  dx run project-G3Yj1vjJ6XG579jbKyjXPGGY:/test/tools/bcftools_isec -f /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/bcftools_isec/bcftools_isec_inputs.json --destination project-G3Yj1vjJ6XG579jbKyjXPGGY:/outputs/