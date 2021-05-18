#!/usr/bin/env python
# coding: utf-8

# Genome Visualizer project

from Bio import SeqIO
import sys

## Function to obtain contig sizes from gbk file, then computes contig locations
## And finally creates a kar file with original contig order.
def main(gbk_filename, output_file):
	
	gbk_file = open(gbk_filename,"r")

	ends = []

	for record in SeqIO.parse(gbk_file, "genbank"):
		location = record.features[0].location
		location = str(location)[1:-4].split(":")
		init = int(location[0])
		end = int(location[1])
		ends.append(end)
	
	# Saving original order
	dic_ends = {i:ends[i] for i in range(len(ends))}
	# Sorting by contig size 
	dic_sorted = {k: v for k, v in sorted(dic_ends.items(), key=lambda item: item[1])}
	
	ends_sort = list(dic_sorted.values())
	new_ends = []
	
	inits = []
	init = 0
	end = ends_sort[0]
	for i in range(len(ends)-1):
		inits.append(init)
		new_ends.append(end)
		init = end + 1
		end = init + ends_sort[i+1]
	inits.append(init)
	new_ends.append(end)
		
	chrx = 1
	chrx = str(chrx).zfill(2)
	lines = []
	
	lines.append("chr - chr"+chrx+" 1 0 "+str(ends_sort[-1])+" black\n")
	for i in range(len(inits)):
		if len(inits) == 1:
			break
		if i%2 == 0:
			color = " white\n"
		else:
			color = " black\n"
		band = str(i+1).zfill(2)
		value = dic_sorted.get(i)
		index = list(dic_sorted.values()).index(value)
		line = "band chr"+chrx+" band"+band+" band"+band+" "+str(inits[index])+" "+str(new_ends[index])+color
		lines.append(line)

	with open(output_file, 'w') as output:
		output.writelines(lines)
		print(output_file+" created succesfully.")

def help():
	print("Run: python create_kar.py -i <input genbank file path> -o <output kar file path>")			

if __name__ == '__main__':
	
	if len(sys.argv) == 1:
		help()
	elif len(sys.argv) == 2 and sys.argv[1] in ['-h', '--h', '-help', '--help', '-H', '--H']:
		help()
	elif len(sys.argv) == 5:
		if sys.argv[1] == '-i' and sys.argv[3] == '-o':
			gbk_file = sys.argv[2]
			output = sys.argv[4]
		elif sys.argv[3] == '-i' and sys.argv[1] == '-o':
			gbk_file = sys.argv[4]
			output = sys.argv[2]
		
		main(gbk_file, output)
	else:
		help()
	
	
