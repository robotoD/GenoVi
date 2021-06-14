#!/usr/bin/env python
# coding: utf-8

# Genome Visualizer project

from Bio import SeqIO
import sys
import numpy as np
import csv

## Function to obtain contig sizes from gbk file, then computes contig locations
## And finally creates a kar file with original contig order.

def ends_sorted(ends):
	dic_ends = {i:ends[i] for i in range(len(ends))}
	dic_sorted = {k: v for k, v in sorted(dic_ends.items(), key=lambda item: item[1])}
	ends_sort = list(dic_sorted.values())
	indexes = list(dic_sorted.keys())
	return ends_sort, indexes

def create_kar(gbk_filename, output_file):
	
	gbk_file = open(gbk_filename,"r")

	ends = []

	for record in SeqIO.parse(gbk_file, "genbank"):
		location = record.features[0].location
		location = str(location)[1:-4].split(":")
		init = int(location[0])
		end = int(location[1])
		ends.append(end)
	
	# Saving original order
	#dic_ends = {i:ends[i] for i in range(len(ends))}
	# Sorting by contig size 
	#dic_sorted = {k: v for k, v in sorted(dic_ends.items(), key=lambda item: item[1])}
	
	#ends_sort = list(dic_sorted.values())
	
	new_ends = []
	
	inits = []
	init = 0
	end = ends[0]
	for i in range(len(ends)-1):
		inits.append(init)
		new_ends.append(end)
		init = end + 1
		end = end + ends[i+1]
	inits.append(init)
	new_ends.append(end)
		
	chrx = 1
	chrx = str(chrx).zfill(2)
	lines = []
	
	lines.append("chr - chr"+chrx+" 1 0 "+str(new_ends[-1])+" black\n")
	for i in range(len(inits)):
		if len(inits) == 1:
			break
		if i%2 == 0:
			color = " white\n"
		else:
			color = " black\n"
		band = str(i+1).zfill(2)
		line = "band chr"+chrx+" band"+band+" band"+band+" "+str(inits[i])+" "+str(new_ends[i])+color
		lines.append(line)

	with open(output_file, 'w') as output:
		output.writelines(lines)
		print(output_file+" created succesfully.")


	return ends

def create_CDS(gbk_filename, cds_p_output, cds_n_output, sizes):
	
	def new_loc(array, sizes_x):
		new_arr = []
		for i, loc_pos in enumerate(array):
			loc_pos = loc_pos.split(":")
			init = int(loc_pos[0][1:])
			end = int(loc_pos[1][:-1])
			new_arr.append([init+sizes_x[i],end+sizes_x[i]])
		return new_arr
		
	def write_lines(locations, sizes_x, names, output_):
		lines = []
		for i in range(len(locations)):
			#line = ["chr-Node_x_length_"+str(sizes_x[i])+"_cov_x"] + list(map(str, locations[i]))
			line = [names[i]] + list(map(str, locations[i]))
			lines.append(line)
		with open(output_, 'w', newline='') as csvfile:
			writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			writer.writerows(lines)	
			print(output_,"created succesfully.")

	
	gbk_file = open(gbk_filename,"r")
	positives = []
	negatives = []
	sizes_p = []
	sizes_n = []
	names_p = []
	names_n = []
	
	for j, record in enumerate(SeqIO.parse(gbk_file, "genbank")):
		aux_p = []
		aux_n = []
		name = record.name
		size_sum = np.sum(np.array(sizes[:j]))
		
		for feature in record.features:
			if feature.type == "CDS":
				location = str(feature.location)[:-3]
				direction = str(feature.location)[-2:-1]
				if direction == '+':
					positives.append(location)
					if j == 0:
						aux_p.append(0)
					else:
						aux_p.append(size_sum)
					names_p.append(name)
				else:
					negatives.append(location)
					if j == 0:
						aux_n.append(0)
					else:
						aux_n.append(size_sum)
					names_n.append(name)
				
				
		sizes_p = sizes_p + aux_p
		sizes_n = sizes_n + aux_n
	
	new_pos = new_loc(positives, sizes_p)
	new_negs = 	new_loc(negatives, sizes_n)
	
	write_lines(new_pos, sizes_p, names_p, cds_p_output)
	write_lines(new_negs, sizes_n, names_n, cds_n_output)
	
	return
				

def help():
	print("Run")
	print("For only creating KAR file: python create_kar.py -i <input genbank file path> -o <output kar file path>")	
	print("For creating KAR file and CDS positive and CDS negative files: python create_kar.py -i <input genbank file path> -o <output kar file path> -cp <output_CDS_positives> -cn <output_CDS_negatives>")			

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
		
		sizes = create_kar(gbk_file, output)
	elif len(sys.argv) == 9:
		if sys.argv[1] == '-i' and sys.argv[3] == '-o' and sys.argv[5] == '-cp' and sys.argv[7] == '-cn':
			gbk_file = sys.argv[2]
			output = sys.argv[4]
			cds_pos = sys.argv[6]
			cds_neg = sys.argv[8]
		elif sys.argv[3] == '-i' and sys.argv[1] == '-o' and sys.argv[5] == '-cp' and sys.argv[7] == '-cn':
			gbk_file = sys.argv[4]
			output = sys.argv[2]
			cds_pos = sys.argv[6]
			cds_neg = sys.argv[8]
		elif sys.argv[1] == '-i' and sys.argv[3] == '-o' and sys.argv[7] == '-cp' and sys.argv[5] == '-cn':
			gbk_file = sys.argv[4]
			output = sys.argv[2]
			cds_pos = sys.argv[8]
			cds_neg = sys.argv[6]
		elif sys.argv[3] == '-i' and sys.argv[1] == '-o' and sys.argv[7] == '-cp' and sys.argv[5] == '-cn':
			gbk_file = sys.argv[4]
			output = sys.argv[2]
			cds_pos = sys.argv[8]
			cds_neg = sys.argv[6]
		
		sizes = create_kar(gbk_file, output)
		create_CDS(gbk_file, cds_pos, cds_neg, sizes)
	else:
		help()
	
	
