#!/usr/bin/env python
# coding: utf-8

# Genome Visualizer project

from Bio import SeqIO
import sys
import numpy as np
import csv
import argparse as ap
import re
import subprocess
import os
from shutil import which
import pandas as pd

def get_args():
    parser = ap.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str, help="Genbank file path", required=True)
    parser.add_argument("-c", "--complete_genome", action='store_true', help="Indicating if it is a complete genome")
    
    pred_args = parser.add_argument_group('Categories prediction arguments')
    pred_args.add_argument("-gc", "--get_categories", action='store_true', help="Indicating if CDS categories must be predicted")
    pred_args.add_argument("-d", "--divided", action='store_true', help="Indicating if CDS categories must be divided in different files")
        
    kar_args = parser.add_argument_group('KAR generation arguments')
    kar_args.add_argument("-o", "--output_file", type=str, help="Output KAR file path. Default: output.kar", default = "output.kar")
    
    cds_args = parser.add_argument_group('CDSs generation arguments')
    cds_args.add_argument("-cp", "--cds_pos", type=str, help="Positive CDS output file", required = False, default = "cds_pos")
    cds_args.add_argument("-cn", "--cds_neg", type=str, help="Negative CDS output file", required = False, default = "cds_neg")
    
    trna_args = parser.add_argument_group('tRNAs generation arguments')
    trna_args.add_argument("-tp", "--trna_pos", type=str, help="Positive tRNA output file", required = False, default = "trna_pos")
    trna_args.add_argument("-tn", "--trna_neg", type=str, help="Negative tRNA output file", required = False, default = "trna_neg")
    
    args = parser.parse_args()
        
    return args.input_file, args.output_file, args.cds_pos, args.cds_neg, args.trna_pos, args.trna_neg, args.complete_genome, args.get_categories, args.divided


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
		
	chrx = '1'
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


def create_kar_complete(gbk_filename, output_file):
	
	gbk_file = open(gbk_filename,"r")

	ends = []

	for record in SeqIO.parse(gbk_file, "genbank"):
		location = record.features[0].location
		location = str(location)[1:-4].split(":")
		init = int(location[0])
		end = int(location[1])
		ends.append(end)
		
	lines = []
	
	for i in range(len(ends)):
		line1 = "chr - chr"+ str(i+1) +" 1 0 "+ str(ends[i]) +" black\n"
		line2 = "band chr"+ str(i+1)+" band01 band01 0 "+str(ends[i])+" white\n"
		lines.append(line1)
		lines.append(line2)

	with open(output_file, 'w') as output:
		output.writelines(lines)
		print(output_file+" created succesfully.")


	return ends
	

def create_kar_plus(gbk_filename, output_file):
	
	gbk_file = open(gbk_filename,"r")

	ends = []

	for record in SeqIO.parse(gbk_file, "genbank"):
		location = record.features[0].location
		location = str(location)[1:-4].split(":")
		init = int(location[0])
		end = int(location[1])
		ends.append(end)
		
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
		
	lines = []
	
	for i in range(len(ends)):
		line1 = "chr - chr"+ str(i+1) +" 1 0 "+ str(ends[i]) +" black\n"
		line2 = "band chr"+ str(i+1)+" band01 band01 "+str(inits[i])+" "+str(new_ends[i])+" white\n"
		lines.append(line1)
		lines.append(line2)

	with open(output_file, 'w') as output:
		output.writelines(lines)
		print(output_file+" created succesfully.")


	return ends	


def new_loc(array, sizes_x):
	new_arr = []
	for i, loc_pos in enumerate(array):
		loc_pos = loc_pos.split(":")
		init = int(loc_pos[0][1:])
		end = int(loc_pos[1][:-1])
		new_arr.append([init+sizes_x[i],end+sizes_x[i]])
	return new_arr
	
def write_lines(locations, output_, chrx, locus, cogs):
	lines = []		
	for i in range(len(locations)):
		#line = ["chr-Node_x_length_"+str(sizes_x[i])+"_cov_x"] + list(map(str, locations[i]))
		if len(locus) == 0:
			line = ["chr"+chrx[i]] + list(map(str, locations[i]))
		elif len(cogs) == 0:
			line = ["chr"+chrx[i]] + list(map(str, locations[i])) + [locus[i]]
		else:
			line = ["chr"+chrx[i]] + list(map(str, locations[i])) + [locus[i]] + [cogs[i]]
		lines.append(line)
	with open(output_, 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer.writerows(lines)	
		print(output_,"created succesfully.")
		
def write_cog_files(locations, output, chrx, locus, cogs):
	
	cogs_df = pd.DataFrame.from_dict(cogs)
	cogs_df.columns = ["category"]
	cogs_df["location"] = locations
	cogs_df[['loc_init','loc_end']] = pd.DataFrame(cogs_df.location.tolist(), index= cogs_df.index)
	cogs_df["chr"] = chrx
	cogs_df["locus"] = locus
	cogs_df["main"] = cogs_df["category"].str[0]
	categories = cogs_df["main"].unique()
	
	for c in categories:
		lines = []	
		subset = cogs_df.loc[cogs_df["main"] == c]
		subset = subset.sort_values(["chr", "loc_init"], ascending=[True, True])
		for index, row in subset.iterrows():
			line = ["chr"+row["chr"]] + list(map(str, row["location"])) + [row["locus"]] + [row["category"]]
			lines.append(line)
		filename = output.split(".")[-2]+"_"+c+"."+output.split(".")[-1]
		with open(filename, 'w', newline='') as csvfile:
			writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			writer.writerows(lines)	
			print(filename, "created succesfully.")		
	


def create_feature(gbk_filename, p_output, n_output, sizes, feat, cogs_dict=None, divided=False):
	
	#chrx = '1'
	
	gbk_file = open(gbk_filename,"r")
	positives = []
	negatives = []
	sizes_p = []
	sizes_n = []
	names_p = []
	names_n = []
	chrx = '1'
	chrms_p = []
	chrms_n = []
	locus_n = []
	locus_p = []
	cogs_p = []
	cogs_n = []
	
	for j, record in enumerate(SeqIO.parse(gbk_file, "genbank")):
		aux_p = []
		aux_n = []
		name = record.name
		size_sum = np.sum(np.array(sizes[:j]))

		for feature in record.features:
			if feature.type == feat: # "CDS" or "tRNA"
				location = str(feature.location)[:-3].replace("<", "").replace(">", "")
				locus_tag = feature.qualifiers.get("locus_tag")[0]
				if("join" in str(feature.location)):
					locationMatch = re.match(r"join\{\[<?>?(\d+):<?>?\d+\]\(.\),\s(?:\[<?>?\d+:<?>?\d+\]\(.\),\s)*\[<?>?\d+:<?>?(\d+)\]\(.\)\}", str(feature.location)) # Sacamos el primer y ultimo numero que encontremos
					location = "[{}:{}]".format(locationMatch.groups()[0], locationMatch.groups()[1])
				direction = str(feature.location)[-2:-1]
				if direction == '+':
					positives.append(location)
					if j == 0:
						aux_p.append(0)
					else:
						aux_p.append(size_sum)
					names_p.append(name)
					chrms_p.append(chrx)
					if feature.type == "CDS":
						locus_p.append(locus_tag)
						if cogs_dict != None:
							cogs_p.append(cogs_dict.get(locus_tag))
				else:
					negatives.append(location)
					if j == 0:
						aux_n.append(0)
					else:
						aux_n.append(size_sum)
					names_n.append(name)
					chrms_n.append(chrx)
					if feature.type == "CDS":
						locus_n.append(locus_tag)
						if cogs_dict != None:
							cogs_n.append(cogs_dict.get(locus_tag))
				
				
		sizes_p = sizes_p + aux_p
		sizes_n = sizes_n + aux_n
	
	new_pos = new_loc(positives, sizes_p)
	new_negs = new_loc(negatives, sizes_n)
	
	if divided:
		write_cog_files(new_pos, p_output, chrms_p, locus_p, cogs_p)
		write_cog_files(new_negs, n_output, chrms_n, locus_n, cogs_n)
	else:
		write_lines(new_pos, p_output, chrms_p, locus_p, cogs_p)
		write_lines(new_negs, n_output, chrms_n, locus_n, cogs_n)
	
	return

def create_feature_complete(gbk_filename, p_output, n_output, sizes, feat, cogs_dict=None, divided=False):
	
	#chrx = '1'
	
	gbk_file = open(gbk_filename,"r")
	positives = []
	negatives = []
	sizes_p = []
	sizes_n = []
	names_p = []
	names_n = []
	chrms_p = []
	chrms_n = []
	locus_p = []
	locus_n = []
	cogs_p = []
	cogs_n = []
	
	for j, record in enumerate(SeqIO.parse(gbk_file, "genbank")):
		aux_p = []
		aux_n = []
		name = record.name
		size_sum = sizes[j]
		chrx = str(j+1)

		for feature in record.features:
			if feature.type == feat: # "CDS" or "tRNA"

				location = str(feature.location)[:-3]
				direction = str(feature.location)[-2:-1]
				locus_tag = feature.qualifiers.get("locus_tag")[0]
				if direction == '+':
					positives.append(location)
					aux_p.append(size_sum)
					names_p.append(name)
					chrms_p.append(chrx)
					if feature.type == "CDS":
						locus_p.append(locus_tag)
						if cogs_dict != None:
							cogs_p.append(cogs_dict.get(locus_tag))
				else:
					negatives.append(location)
					aux_n.append(size_sum)
					names_n.append(name)
					chrms_n.append(chrx)
					if feature.type == "CDS":
						locus_n.append(locus_tag)
						if cogs_dict != None:
							cogs_n.append(cogs_dict.get(locus_tag))
				
		sizes_p = sizes_p + aux_p
		sizes_n = sizes_n + aux_n
	
	new_pos = new_loc(positives, sizes_p)
	new_negs = 	new_loc(negatives, sizes_n)
	
	if divided:
		write_cog_files(new_pos, p_output, chrms_p, locus_p, cogs_p)
		write_cog_files(new_negs, n_output, chrms_n, locus_n, cogs_n)
	else:
		write_lines(new_pos, p_output, chrms_p, locus_p, cogs_p)
		write_lines(new_negs, n_output, chrms_n, locus_n, cogs_n)
	
	return


def get_categories(gbk_file):
	
	
	# Check if deepnog is installed
	
	if which("deepnog") == None:
		print("Error: Deepnog tool not found. Please follow installation guide on https://github.com/univieCUBE/deepnog")
		return
	
	# Transform gbk into faa
	
	path =  os.getcwd()
	
	try:
		faa_name = gbk_file.split("/")[-1]
	except:
		pass
	faa_name = faa_name.split(".")[0]
	
	output_faa = path + "/" + faa_name + ".faa"
	command1 = "python " + path + "/genbank2faa.py --outputFile " + output_faa + " " + gbk_file
	
	try: 
		print(command1)
		process = subprocess.Popen(command1.split(), stdout=subprocess.PIPE)
		output, error = process.communicate()
		print()
		print("GBK file transformed into faa succesfully. File saved as", output_faa)
	except:
		print("Error when transforming gbk to faa.")
		return
	
	# Predict from deepnog
	
	output_pred = path + "/prediction_deepnog_" + faa_name + ".csv"
	command2 = "deepnog infer " + output_faa + " --out " + output_pred + " -db cog2020 -t 1"
	
	try: 
		print()
		print("Deepnog prediction started")
		print()
		print(command2)
		process = subprocess.Popen(command2.split(), stdout=subprocess.PIPE)
		output, error = process.communicate()
		print()
		print("Deepnog prediction finished succesfully. Predictions saved as", output_pred)
		print()
	except:
		print("Error when predicting categories from " + output_faa + " with deepnog.")
		return
	
	# Cross files 
	
	cogs_df = pd.read_csv(output_pred, sep=',', usecols=['sequence_id', 'prediction'])
	tab_file = path + "/dataset/cog-20.def.tab"
	cogs_df.columns = ['id', 'cog']
	tab_df = pd.read_csv(tab_file, header=None, sep='\t', usecols=[0,1])
	tab_df.columns = ['cog', 'category']
	merge_df = pd.merge(cogs_df, tab_df, on=["cog"])
	
	cogs_dict = merge_df.set_index('id')['category'].to_dict()
	
	return cogs_dict
	
		
if __name__ == '__main__':
	
	gbk_file, output, cds_pos, cds_neg, trna_pos, trna_neg, com_gen, get_cats, divided = get_args()[:]
	
	flag = True
	
	if cds_pos == "cds_pos" and cds_neg == "cds_neg" and get_cats == True:
		print("Error: Categories can only be predicted for CDS. Please enter output file paths for both CDS.") 
		flag = False
	elif get_cats == True:
		cogs_dict = get_categories(gbk_file)
	else:
		cogs_dict = None
		
	if divided == True and get_cats == False:
		print("Error: Division of categories is part of the categories prediction process. Include --get_categories in your command.") 
		flag = False
	
	if flag == True:
		if cds_pos == "cds_pos" and cds_neg == "cds_neg" and trna_pos == "trna_pos" and trna_neg == "trna_neg" and com_gen == False:
			sizes = create_kar(gbk_file, output)
		elif cds_pos == "cds_pos" and cds_neg == "cds_neg" and trna_pos == "trna_pos" and trna_neg == "trna_neg" and com_gen == True:
			sizes = create_kar_plus(gbk_file, output)

		elif (cds_pos == "cds_pos" and cds_neg != "cds_neg") or (cds_neg == "cds_neg" and cds_pos != "cds_pos"):
			print("Error: Please enter an output file path for both CDS positives and CDS negatives.") 
		elif (trna_pos == "trna_pos" and trna_neg != "trna_neg") or (trna_neg == "trna_neg" and trna_pos != "trna_pos"):
			print("Error: Please enter an output file path for both tRNA positives and tRNA negatives.") 
			
		elif cds_pos == "cds_pos" and cds_neg == "cds_neg" and trna_pos != "trna_pos" and trna_neg != "trna_neg" and com_gen == True:
			sizes = create_kar_plus(gbk_file, output)
			create_feature_complete(gbk_file, trna_pos, trna_neg, sizes, "tRNA")
		elif cds_pos == "cds_pos" and cds_neg == "cds_neg" and trna_pos != "trna_pos" and trna_neg != "trna_neg" and com_gen == False:
			sizes = create_kar(gbk_file, output)
			create_feature(gbk_file, trna_pos, trna_neg, sizes, "tRNA")
			
		elif cds_pos != "cds_pos" and cds_neg != "cds_neg" and trna_pos == "trna_pos" and trna_neg == "trna_neg" and com_gen == True:
			sizes = create_kar_plus(gbk_file, output)
			create_feature_complete(gbk_file, cds_pos, cds_neg, sizes, "CDS", cogs_dict, divided)
		elif cds_pos != "cds_pos" and cds_neg != "cds_neg" and trna_pos == "trna_pos" and trna_neg == "trna_neg" and com_gen == False:
			sizes = create_kar(gbk_file, output)
			create_feature(gbk_file, cds_pos, cds_neg, sizes, "CDS", cogs_dict, divided)
			
		elif cds_pos != "cds_pos" and cds_neg != "cds_neg" and trna_pos != "trna_pos" and trna_neg != "trna_neg" and com_gen == True:
			sizes = create_kar_plus(gbk_file, output)
			create_feature_complete(gbk_file, cds_pos, cds_neg, sizes, "CDS", cogs_dict, divided)
			create_feature_complete(gbk_file, trna_pos, trna_neg, sizes, "tRNA")
		else:	
			sizes = create_kar(gbk_file, output)
			create_feature(gbk_file, cds_pos, cds_neg, sizes, "CDS", cogs_dict, divided)
			create_feature(gbk_file, trna_pos, trna_neg, sizes, "tRNA")

	
	
