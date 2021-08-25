#! /usr/bin/env python3


import json
import argparse

parser = argparse.ArgumentParser(description='CRAM entry.')
parser.add_argument('eid', metavar='N', type=str,
                    help='key')
parser.add_argument('cram', metavar='N', type=str,
                    help='key')
parser.add_argument('cram_index', metavar='N', type=str,
                    help='key')
parser.add_argument('folder', metavar='N', type=str,
                    help='output folder with template to name new input file')
parser.add_argument('folder_prefix', metavar='N', type=str,
                    help='output folder with template to name new input file')
# args = parser.parse_args(['1036101','"file-Fyf3kj8J7YXfQX8B3kyz1YYY"',
#                           "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/Mutect_Vardict_CH_NO_BQSR"])
args = parser.parse_args()
args.cram = args.cram.strip("\"")
args.cram_index = args.cram_index.strip("\"")
with open(args.folder + "/input_json/template2.json", "r") as f:
    my_json = json.load(f)
    
for k, v in my_json.items():
    if k.endswith("bam"):
        my_json[k]['$dnanexus_link']['id'] = args.cram
    elif k.endswith("bam_index"):
        my_json[k]['$dnanexus_link']['id'] = args.cram_index


with open(args.folder + "/input_json/{}/{}.json".format(args.folder_prefix, args.eid), 'w') as outfile:
    json.dump(my_json, outfile, indent=4)