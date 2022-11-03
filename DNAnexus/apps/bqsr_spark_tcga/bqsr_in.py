#! /usr/bin/env python3


import json
import argparse

parser = argparse.ArgumentParser(description='CRAM entry.')
parser.add_argument('eid', metavar='N', type=str,
                    help='key')
parser.add_argument('bam', metavar='N', type=str,
                    help='key')
parser.add_argument('bam_index', metavar='N', type=str,
                    help='key')
parser.add_argument('app_folder', metavar='N', type=str,
                    help='output folder with template to name new input file')


args = parser.parse_args()
args.bam = args.bam.strip("\"")
args.bam_index = args.bam_index.strip("\"")
with open(args.app_folder + "/bqsr_spark_template.json", "r") as f:
    my_json = json.load(f)
    
for k, v in my_json.items():
    if k.endswith("bam"):
        my_json[k]['$dnanexus_link']['id'] = args.bam
    elif k.endswith("bam_index"):
        my_json[k]['$dnanexus_link']['id'] = args.bam_index


with open(args.app_folder + "/{}.json".format(args.eid), 'w') as outfile:
    json.dump(my_json, outfile, indent=4)