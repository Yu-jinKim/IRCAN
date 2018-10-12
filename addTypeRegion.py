#!/usr/bin/env python3

usage="\
addTypeRegion.py -m FILE ...\
	\n\n\tAdd type of capture regions to the CapStarr-Seq master file.\n\
"
version="addTypeRegion.py version 1.0"

import sys
import argparse
import fileinput

from signal import signal, SIGPIPE, SIG_DFL

signal(SIGPIPE, SIG_DFL) 

parser = argparse.ArgumentParser()
parser.add_argument("-m","--master", type=argparse.FileType('r'))
args = parser.parse_args()

if not args.master :
	sys.exit("No master file given, please use addTypeRegion.py -h")
elif args.master=="stdin" or args.master=="-":
	fi = fileinput.FileInput(args.master)
else:
	fi = args.master

hardData={}

for index,line in enumerate(fi):
	line = line.split()
	location=tuple(line[0:8])
	if location not in hardData:
		hardData[location]=["_".join(line[8:])]
	else:
		if line[11] not in hardData[location]:
			hardData[location].append("_".join(line[8:]))

for key in hardData:
	print("\t".join(key)+"\t"+"/".join(hardData[key]))