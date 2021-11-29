#! /usr/bin/env python3


import json
import argparse

parser = argparse.ArgumentParser(description='CRAM entry.')
parser.add_argument('name', metavar='N', type=str,
                    help='key')
parser.add_argument('tumor_bam', metavar='N', type=str,
                    help='key')
parser.add_argument('tumor_bai', metavar='N', type=str,
                    help='key')
parser.add_argument('normal_bam', metavar='N', type=str,
                    help='key')
parser.add_argument('normal_bai', metavar='N', type=str,
                    help='key')
parser.add_argument('folder', metavar='N', type=str,
                    help='key')


# =============================================================================
# /Users/brian/Bolton/UKBB/docs/DNAnexus/apps/archer_in.py 218281_1_10.json "file-G5kbfj0JQ288j3YZ0K4YYg7F" "file-G5kbfkQJQ28Jg8Py0Gyfx78K" "file-G5kf6G0JQ28JJbQjPg42YVXB" "file-G5kf6JQJQ28BqGByJ63Yxz1X" 
# 
# args = parser.parse_args(['218281_1_10.json','"file-G5kbfj0JQ288j3YZ0K4YYg7F"','"file-G5kbfkQJQ28Jg8Py0Gyfx78K"','"file-G5kf6G0JQ28JJbQjPg42YVXB"','"file-G5kf6JQJQ28BqGByJ63Yxz1X"',
#                           "/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/Archer"])
# =============================================================================
args = parser.parse_args()
args.tumor_bam = args.tumor_bam.strip("\"")
args.tumor_bai = args.tumor_bai.strip("\"")
args.normal_bam = args.normal_bam.strip("\"")
args.normal_bai = args.normal_bai.strip("\"")
with open("/Users/brian/Bolton/UKBB/docs/DNAnexus/apps/workflows/FINAL/batch/Archer/archer_template.json", "r") as f:
    my_json = json.load(f)
    
for k, v in my_json.items():
    if k.split(".")[1] == "bam":
        my_json[k]['$dnanexus_link']['id'] = args.tumor_bam
    elif k.split(".")[1] == "bam_index" :
        my_json[k]['$dnanexus_link']['id'] = args.tumor_bai
    elif k.split(".")[1] == "normal_bam" :
        my_json[k]['$dnanexus_link']['id'] = args.normal_bam
    elif k.split(".")[1] == "normal_bai" :
        my_json[k]['$dnanexus_link']['id'] = args.normal_bai


with open(args.folder + "/" + args.name, 'w') as outfile:
    json.dump(my_json, outfile, indent=4)