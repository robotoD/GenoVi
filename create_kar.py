#!/usr/bin/env python
# coding: utf-8

# Genome Visualizer project

from Bio import SeqIO
import sys

def main(gbk_filename, output_file):
	
	gbk_file = open(gbk_filename,"r")

	ends = []

	for record in SeqIO.parse(gbk_file, "genbank"):
		location = record.features[0].location
		location = str(location)[1:-4].split(":")
		init = int(location[0])
		end = int(location[1])
		ends.append(end)

	ends.sort()

	inits = []
	init = 0
	for end in ends:
		inits.append(init)
		init = end + 1

	chrx = 1
	chrx = str(chrx).zfill(2)
	lines = []

	lines.append("chr - chr"+chrx+" 1 0 "+str(ends[-1])+" black\n")
	for i in range(len(inits)):
		if len(inits) == 1:
			break
		if i%2 == 0:
			color = " white\n"
		else:
			color = " black\n"
		band = str(i+1).zfill(2)
		line = "band chr"+chrx+" band"+band+" band"+band+" "+str(inits[i])+" "+str(ends[i])+color
		lines.append(line)

	with open(output_file, 'w') as output:
		output.writelines(lines)
		

if __name__ == '__main__':
	gbk_file = sys.argv[1]
	output = sys.argv[2]
	main(gbk_file, output)
	
	
