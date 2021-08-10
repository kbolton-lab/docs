#! /usr/bin/env python3

import sys
import json
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('key', metavar='N', type=str,
                    help='key')
parser.add_argument('indent', metavar='N', type=int,
                    help='indent')
args = parser.parse_args()
# for line in sys.stdin:
#     sys.stdout.write(line)
txt = [{"$dnanexus_link": line.strip()} for line in sys.stdin]


# with open('file_id.txt') as f:
#     txt = [{'$dnanexus_link': line.strip()} for line in f]

#print(txt)

d = {args.key: txt}
print(json.dumps(d, indent=args.indent))