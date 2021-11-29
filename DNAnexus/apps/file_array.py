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

txt = [{"$dnanexus_link": line.strip()} for line in sys.stdin]

d = {args.key: txt}
print(json.dumps(d, indent=args.indent))