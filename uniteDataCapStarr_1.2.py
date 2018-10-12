#!/usr/bin/env python3

usage="\
uniteDataCapStarr.py -m MASTER FILE -f FILE FILE ...\
	\n\n\tFrom CapStarrSeq data, outputs the data for each location \
	\n\tChecks values of nseqs per clone and gets the highest one but gives it to every location of the clone \
	\n\tInserts NAs for files that didn't contribute to the location \
"
version="uniteDataCapStarr.py version 1.1"

import sys
import argparse
import gzip
from signal import signal, SIGPIPE, SIG_DFL
from collections import Counter

def chrSort(keys):
	"""
		Takes an unordered list of keys (location from bed file)
		Returns ordered list according to position
	"""
	chrNb=[]
	chrStr=[]
	for key in keys:															# for every key 
		key=key.split()															# split the line for latter manipulation
		try:
			chrNb.append([int(key[0]),int(key[1]),int(key[2]),key[3]])			# if the first element is an int: becomes an int and append the location 
		except ValueError:
			chrStr.append([str(key[0]),int(key[1]),int(key[2]),key[3]])			# if not: becomes a string and append the location
	chrNb.sort(key=lambda x: (x[0], x[1], x[2], x[3]))							# ordering of chrNb list
	chrStr.sort(key=lambda x: (x[0], x[1], x[2], x[3]))							# ordering of chrStr list 
	return chrNb+chrStr															# returns the 2 lists concatenated

def flatten(x, sep=" "):
	"""
		Returns flatten list (list of lists --> list of str)
	"""
	res=[]
	for ele in x:
		if (ele.__class__ is list):
			ele = sep.join(str(x) for x in ele)
		res.append(ele)
	return res

signal(SIGPIPE, SIG_DFL) 

parser = argparse.ArgumentParser()

parser.add_argument("-m","--master", type=argparse.FileType('r'))
parser.add_argument("-f","--files", type=argparse.FileType('r'), nargs='+')
parser.add_argument("--version", action='store_true')
args = parser.parse_args()

if args.version:
	sys.exit(version)

if args.files:
	if len(args.files) < 1:
		sys.exit("More than 1 file needed\n")
else:
	sys.exit("No files given, please use uniteDataCapStarr_1.1.py -h")

if not args.master :
	sys.exit("No master file given, please use uniteDataCapStarr_1.1.py -h")
elif args.master=="stdin" or args.master=="-":
	fi = fileinput.FileInput(args.master)
else:
	fi = args.master

loc_data={}
counterFile={}
loc_clone={}
clone_data={}
nbFiles=len(args.files)
i=0
j=0

for line in fi:																# for every line in the master file
	line=line.split()
	if line[0]=="#" or line[0]=="Track":									# ignore track and comment lines
		pass
	else:
		d=[]
		location=line[0]+" "+line[1]+" "+line[2]+" "+line[5]				# get location (BED files)
		for x in range(6,len(line)):										# get data (nb seqs, nb locs..)
			d.append(line[x])
		data=" ".join(d)
		if location not in loc_clone:										# new location in loc_clone
			counterFile[location] = []											# initiate dict at location for future annotation of origin of data
			loc_clone[location]=line[3]											# location -> clone
			loc_data[location]={0:data}											# location -> file -> data

fi.close()

for file in args.files:														# processing the files
	dataset_clone=[]
	i+=1
	line = file.readline()
	while line:																# for every line in the file
		line=line.split()
		if line[0]=="#" or line[0]=="Track":								# ignore track and comment lines
			pass
		else:
			d=[]
			location=line[0]+" "+line[1]+" "+line[2]+" "+line[5]			# get location (BED files)
			for x in range(6,len(line),2):									# get data (nb seqs, nb locs..)
				d.append(int(line[x]))
			dataset_clone.append(line[3])
			if location in loc_clone:
				# print("\t".join(line))
				if loc_clone[location] in clone_data and i in clone_data[loc_clone[location]] and line[3] not in dataset_clone:				# clone already exists and in the same file
					clone_data[loc_clone[location]][i]+=int(d[0])												# clone -> file -> stored data + data[0]
				elif loc_clone[location] in clone_data and i not in clone_data[loc_clone[location]]:		# clone already exists and different file
					clone_data[loc_clone[location]][i]=d[0]														# clone -> file -> data[0]
				else:																						# first time encountering the clone
					clone_data[loc_clone[location]] = {i:d[0]}													# clone -> file -> data[0]
				loc_data[location][i]=d																		# location -> file -> data
		# print("\t".join(line))
		line = file.readline()
	file.close()

lengthData=len(d)

insertNA=lengthData*"NA "
insertNA=insertNA.strip()

for k in loc_data:															# insert NA for non contributing files for every location
	j=1
	while j!=nbFiles+1:
		if loc_clone[k] in clone_data and j not in loc_data[k]:				# if clone at this location exists and file didn't contribute
			loc_data[k][j]=["NA","NA"]											# insertion
		if loc_clone[k] in clone_data and j in clone_data[loc_clone[k]]:	# if clone at this location exists and file contributed towards the clone -> file -> data
			loc_data[k][j][0]=clone_data[loc_clone[k]][j]						# remplace value to the sum of the clone
		j+=1

keys = list(loc_clone)														# get the keys (locations) from the location => clone name dictionary
orderedKeys = chrSort(keys)													# get the ordered keys

for key in orderedKeys:
	i=0
	data=[]
	location=" ".join(str(ele) for ele in key).split()
	key=" ".join(str(ele) for ele in key)
	while i <= nbFiles:
		data.append(loc_data[key][i])
		i+=1
	data="\t".join(" ".join(flatten(data)).split())
	# print(data)
	print(str(location[0].strip()+"\t"+location[1].strip()+"\t"+location[2].strip()+"\t"+"".join(loc_clone[key]).strip()+"\t1\t"+location[3].strip()+"\t"+data))
