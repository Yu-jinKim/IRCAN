#!/usr/bin/env python3

"""
uniqCloneName version 1.0:

	If you use a stdin, use the "-" or "stdin"
	
	Takes one bed file argument (6th column NEED to be the nb of aligning locations in capRegions)

	This script will assign a unique clone name to every seqID that align to a location
	Computes nb of unique seqID for each clone name (6h column)
	Computes nb of unique locations for each clone name (7th column)
	Computes nb of alignments in capReg for the location (8th column)

	Outputs sorted bed file for every unique location

	Example output:
	chr chrStart chrEnd clone_name placeholderScore strand nbUniqSeqID_perCloneName nbUniqLocations_perCloneName nbAlignmentCapReg
"""

import sys, argparse, gzip, fileinput, os
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL) 

def chrSort(keys):
	"""
		Takes an unordered list of keys (location from bed file)
		Returns ordered list
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

seq_clone={} 		# seqID => clone name
loc_clone={} 		# location => clone name
countSeq={} 		# clone name => nb uniq seqID
countLoc={} 		# clone name => nb loc
countCapReg={}		# clone name => nb align capReg
seq_loc={}
loc_seq={}
i=0					# for clone name

parser = argparse.ArgumentParser(description="""Outputs sorted BED file with unique clone name,
												 nb uniq SeqID/nb uniq locations per clone as well as nb 
												aligning locations in capture regions per clone""")
parser.add_argument("file", help="bed_file_input")
args = parser.parse_args()

if args.file.lower().endswith('.gz'):											# handling gzip files
	data=gzip.open(args.file,'r')
	encodedData=data.readlines()
	data.close()
	fi=[]
	for ele in encodedData:
		fi.append(ele.decode())
else:																			# handling normal files and stdin
	if args.file=="stdin":
		args.file="-"
	fi = fileinput.FileInput(args.file)

for line in fi:
	line=line.strip("\n").split("\t")												# 	formatting for processing
	seq_id=line[3]																	# 	get the seq ID
	location=[line[0][3:],line[1],line[2],line[5]]									# 	get the location (chr,chrStart,chrEnd,strand)
	location=" ".join(location)														# 	convert in string for use as key in dict
	capReg=int(line[6])																#	get the nb alignment in capReg for the line
	if seq_id not in seq_clone.keys() and location not in loc_clone.keys():			# 	if both seq_ID and location have not been encountered yet :
		i+=1																		# 		increment i (for the name of the clone)
		clone_name="clone"+str(i)													# 		create a new clone name
		seq_clone[seq_id]=clone_name												# 		seqID -> clone name
		loc_clone[location]=clone_name												# 		location -> clone name
		countLoc[clone_name]=1														# 		start counter for nb locations
		countSeq[clone_name]=1														# 		start counter for nb seqID
		countCapReg[clone_name]=capReg												#		initiate the dict with the nb of alignment to capReg
		seq_loc[seq_id]=[location]
		loc_seq[location]=[seq_id]
	elif seq_id not in seq_clone.keys() and location in loc_clone.keys():			# 	if the seqID has not been encountered but location has:
		seq_clone[seq_id]=loc_clone[location]										# 		assign clone name of the location to the seqID encountered
		clone_name=loc_clone[location]												# 		change clone name (for counting)
		countSeq[clone_name]+=1														# 		increment nb uniq seqID
		seq_loc[seq_id]=[]
		loc_seq[location].append(seq_id)
	elif seq_id in seq_clone.keys() and location not in loc_clone.keys():			# 	if the seqID has been encountered but not the location:
		loc_clone[location]=seq_clone[seq_id]										# 		assign the clone name of the seqID to the location
		clone_name=seq_clone[seq_id]												# 		change clone name (for counting)
		countLoc[clone_name]+=1														# 		increment nb uniq location
		countCapReg[clone_name]+=capReg												#		sum up the nb alignment capReg
		seq_loc[seq_id].append(location)
		loc_seq[location]=[]
	elif seq_id in seq_clone.keys() and location in loc_clone.keys():							# if both seqID and location have been encountered:
		oldCloneName=seq_clone[seq_id]															# 	save the old clone name
		newCloneName=loc_clone[location]														# 	save the new clone name
		if seq_clone[seq_id]!=loc_clone[location]:												# 	if clone name do not match: --> ambiguity on clone name
			countSeq[newCloneName]+=countSeq[oldCloneName]										# 		add the seq ID of the old clone name to the new one
			countLoc[newCloneName]+=countLoc[oldCloneName]										# 		add the locations of the old clone name to the new one
			for loc in seq_loc[seq_id]:
				loc_clone[loc]=newCloneName
				for seq in loc_seq[loc]:
					seq_clone[seq]=newCloneName

keys = list(loc_clone.keys())														# get the keys (locations) from the location => clone name dictionary
orderedKeys = chrSort(keys)															# get the ordered keys

print("Track name=\""+args.file+"\" description=\"created with uniqCloneName_1.3.py\"")
for key in orderedKeys:																		# printing
	position=" ".join(str(ele) for ele in key)
	print("chr"+str(key[0]),"\t",															# chrNb
			str(key[1]),"\t",																# chrStart
			str(key[2]),"\t",																# chrEnd
			loc_clone[position],"\t",														# clone name
			"1","\t",																		# placeholder Score
			key[3],"\t",																	# strand
			countSeq[loc_clone[position]],"\t",												# nb of uniq SeqID
			countLoc[loc_clone[position]],"\t",												# nb of uniq location
			countCapReg[loc_clone[position]])												# nb of alignment capReg for this location
