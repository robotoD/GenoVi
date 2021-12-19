#!/usr/bin/env python
# coding: utf-8

# GenoVi is a pipeline that generates circular maps for bacterial (complete or non-complete)
# genomes using Circos software. It also allows the user to annotate COG classifications
# through DeepNOG predictions.
# 
# GenoVi is under a BY-NC-SA Creative Commons License, Please cite. Cumsille et al., 2021
# You may remix, tweak, and build upon this work even for commercial purposes, as long as
# you credit this work and license your new creations under the identical terms.
# 
# Developed by Andres Cumsille, Andrea Rodriguez, Roberto E. Duran & Vicente Saona Urmeneta
# For any code related query, contact: andrea.rodriguezdelherbe@rdm.ox.ac.uk, vicente.saona@sansano.usm.cl.


from Bio import SeqIO
import numpy as np
import csv
import argparse as ap
import re
import subprocess
import os
from shutil import which
import pandas as pd
import scripts.genbank2faa as genbank2faa
global seek


# Parse user arguments
def get_args():
    parser = ap.ArgumentParser()
    parser.add_argument("input_file", help="Genbank file path")
    
    parser.add_argument("-o", "--output_folder", type=str, help="Output folder path. By default it will take the name of the gbk file", required = False, default="")
    
    #cds_args = parser.add_argument_group('CDSs generation arguments')
    #cds_args.add_argument("-cp", "--cds_pos", type=str, help="Positive CDS output file", required = False, default = "cds_pos")
    #cds_args.add_argument("-cn", "--cds_neg", type=str, help="Negative CDS output file", required = False, default = "cds_neg")
    parser.add_argument("-cds", "--cds", action='store_true', help="CDS files (positive and negative) will be created", required = False)
    
    #trna_args = parser.add_argument_group('tRNAs generation arguments')
    #trna_args.add_argument("-tp", "--trna_pos", type=str, help="Positive tRNA output file", required = False, default = "trna_pos")
    #trna_args.add_argument("-tn", "--trna_neg", type=str, help="Negative tRNA output file", required = False, default = "trna_neg")
    parser.add_argument("-trna", "--trna", action='store_true', help="tRNA files (positive and negative) will be created", required = False)
    parser.add_argument("-rrna", "--rrna", action='store_true', help="rRNA files (positive and negative) will be created", required = False)    
        
    pred_args = parser.add_argument_group('Categories prediction arguments')
    pred_args.add_argument("-gc", "--get_categories", action='store_true', help="Indicating if CDS categories must be predicted")
    pred_args.add_argument("-d", "--divided", action='store_true', help="Indicating if CDS categories must be divided in different files")
    
    parser.add_argument("-c", "--complete_genome", action='store_true', help="Indicating if it is a complete genome")
    
    args = parser.parse_args()
        
    return args.input_file, args.output_folder, args.cds, args.trna, args.rrna, args.get_categories, args.divided, args.complete_genome

# Recursive function that looks for genovi folder
def listdir_r(dirpath, folder, seek):
	for path in os.listdir(dirpath):
		rpath = os.path.join(dirpath, path)
		if "genovi" in rpath:
			if os.path.isdir(rpath):
				subdirs = listdir_r(rpath, folder, seek)
				
			else:
				split = rpath.split('/')
				if split[-1] == "genovi":
					direct = '/'.join(split[:-1])
					if os.path.isdir(direct + "/input_test"):
						seek.append(direct)
		else:
			if os.path.isdir(rpath):
				subdirs = listdir_r(rpath, folder, seek)		
		
	return seek
	

# Function to obtain contig sizes from gbk file, then computes contig locations
# And finally creates a kar file with original contig order.
def ends_sorted(ends):
	dic_ends = {i:ends[i] for i in range(len(ends))}
	dic_sorted = {k: v for k, v in sorted(dic_ends.items(), key=lambda item: item[1])}
	ends_sort = list(dic_sorted.values())
	indexes = list(dic_sorted.keys())
	return ends_sort, indexes

# Fuction for creating base KAR file that defines contig bands.
# It considers that the genome is not complete.
def create_kar(gbk_filename, output_folder, complete):
	
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
	
	if complete == False:
		for i in range(len(ends)):
			line1 = "chr - chr"+ str(i+1) +" 1 " +str(inits[i])+" "+str(new_ends[i])+" black\n"
			if i%2 == 0:
				color = " white\n"
			else:
				color = " black\n"
			line2 = "band chr"+ str(i+1)+" band01 band01 "+str(inits[i])+" "+str(new_ends[i])+color
			lines.append(line1)
			lines.append(line2)

		output_file = output_folder + "_bands.kar"
		with open(output_file, 'w') as output:
			output.writelines(lines)
			print(output_file+" created succesfully.")


	return ends, inits, new_ends		

# Fuction for creating base KAR file that defines contig bands.
# It considers that the genome is complete.
def create_kar_complete(output_folder, k, init, end):
		
	lines = []
	
	line1 = "chr - chr"+ str(k+1) +" 1 " +str(init)+" "+str(end)+" black\n"
	line2 = "band chr"+ str(k+1)+" band01 band01 "+str(init)+" "+str(end)+" white\n"
	lines.append(line1)
	lines.append(line2)

	output_file = output_folder + "_bands.kar"
	with open(output_file, 'w') as output:
		output.writelines(lines)
		print(output_file+" created succesfully.")


	return

# Calculating locations of contigs.
def new_loc(array, sizes_x):
	new_arr = []
	for i, loc_pos in enumerate(array):
		loc_pos = loc_pos.split(":")
		init = int(loc_pos[0][1:])
		end = int(loc_pos[1][:-1])
		new_arr.append([init+sizes_x[i],end+sizes_x[i]])
		#new_arr.append([init,end])
	return new_arr
	
# Writting locations to file.
def write_lines(locations, output_, chrx, locus, cogs, verbose = False):
	lines = []
	#print(locations)		
	for i in range(len(locations)):
		#line = ["chr-Node_x_length_"+str(sizes_x[i])+"_cov_x"] + list(map(str, locations[i]))
		if len(locus) == 0:
			line = ["chr"+chrx[i]] + list(map(str, locations[i]))
		elif len(cogs) == 0:
			line = ["chr"+chrx[i]] + list(map(str, locations[i]))# + [locus[i]]
		else:
			line = ["chr"+chrx[i]] + list(map(str, locations[i]))# + [locus[i]] + [cogs[i]]
		lines.append(line)
	with open(output_, 'w', newline='') as csvfile:
		writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
		writer.writerows(lines)
		if verbose:
			print(output_,"created succesfully.")

def write_cog_files(locations, output, chrx, locus, cogs, verbose = False):
	
	if len(cogs) == 0:
		return
	cogs_df = pd.DataFrame.from_dict(cogs)
	cogs_df.columns = ["category"]
	cogs_df["location"] = locations
	cogs_df[['loc_init','loc_end']] = pd.DataFrame(cogs_df.location.tolist(), index= cogs_df.index)
	cogs_df["chr"] = chrx
	cogs_df["locus"] = locus
	cogs_df["main"] = cogs_df["category"].str[0]
	categories = map(str, cogs_df["main"].unique())
	
	for c in categories:
		lines = []
		subset = cogs_df.loc[cogs_df["main"] == c]
		subset = subset.sort_values(["chr", "loc_init"], ascending=[True, True])
		for index, row in subset.iterrows():
			line = ["chr"+row["chr"]] + list(map(str, row["location"]))
			lines.append(line)
		filename = output.split(".")[-2]+"_"+c+"."+output.split(".")[-1]
		with open(filename, 'w', newline='') as csvfile:
			writer = csv.writer(csvfile, delimiter='\t', quotechar='|', quoting=csv.QUOTE_MINIMAL)
			writer.writerows(lines)
			if verbose:
				print(filename, "created succesfully.")
	

# Creates feature (CDS, tRNAm, or rRNA) files for CIRCOS.
# It considers that the genome is not complete.
def create_feature(gbk_filename, output, sizes, feat, cogs_dict=None, divided=False, verbose = False):
	
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
		size_sum = np.sum(np.array(sizes[:j]))
		chrx = str(j+1)

		for feature in record.features:
			if feature.type == feat: # "CDS" or "tRNA"

				location = str(feature.location)[:-3].replace("<", "").replace(">", "")
				if("join" in str(feature.location)):
					locationMatch = re.match(r"join\{\[<?>?(\d+):<?>?(\d+)\]\(.\),\s(?:\[<?>?\d+:<?>?\d+\]\(.\),\s)*\[<?>?\d+:<?>?(\d+)\]\(.\)\}", str(feature.location)) # Sacamos el primer y ultimo numero que encontremos
					if 0 in feature.location:
						location = "[{}:{}]".format(locationMatch.groups()[0], locationMatch.groups()[1])
					else:
						location = "[{}:{}]".format(locationMatch.groups()[0], locationMatch.groups()[2])
				direction = str(feature.location)[-2:-1]
				locus_tag = feature.qualifiers.get("locus_tag")[0]
				if direction == '+':
					if j == 0:
						aux_p.append(0)
					else:
						aux_p.append(size_sum)
					positives.append(location)
					names_p.append(name)
					chrms_p.append(chrx)
					if feature.type == "CDS":
						locus_p.append(locus_tag)
						if cogs_dict != None:
							cogs_p.append(cogs_dict.get(locus_tag))
				else:
					if j == 0:
						aux_n.append(0)
					else:
						aux_n.append(size_sum)
					negatives.append(location)
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
	
	p_output = output + "_" + feat + "_pos.txt"
	n_output = output + "_" + feat + "_neg.txt"
	
	if divided:
		write_cog_files(new_pos, p_output, chrms_p, locus_p, cogs_p, verbose = verbose)
		write_cog_files(new_negs, n_output, chrms_n, locus_n, cogs_n, verbose = verbose)
	else:
		write_lines(new_pos, p_output, chrms_p, locus_p, cogs_p, verbose = verbose)
		write_lines(new_negs, n_output, chrms_n, locus_n, cogs_n, verbose = verbose)
	
	return(cogs_p, cogs_n)
	

# Creates feature (CDS, tRNAm, or rRNA) files for CIRCOS.
# It considers that the genome is complete.
def create_feature_complete(gbk_filename, output, sizes, j, feat, cogs_dict=None, divided=False):
	
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
	size_sum = np.sum(np.array(sizes[:j]))
	
	for record in SeqIO.parse(gbk_file, "genbank"):
		aux_p = []
		aux_n = []
		name = record.name
		chrx = str(j+1)

		for feature in record.features:
			if feature.type == feat: # "CDS" or "tRNA"

				location = str(feature.location)[:-3].replace("<", "").replace(">", "")
				if("join" in str(feature.location)):
					locationMatch = re.match(r"join\{\[<?>?(\d+):<?>?\d+\]\(.\),\s(?:\[<?>?\d+:<?>?\d+\]\(.\),\s)*\[<?>?\d+:<?>?(\d+)\]\(.\)\}", str(feature.location)) # Sacamos el primer y ultimo numero que encontremos
					location = "[{}:{}]".format(locationMatch.groups()[0], locationMatch.groups()[1])
				direction = str(feature.location)[-2:-1]
				locus_tag = feature.qualifiers.get("locus_tag")[0]
				if direction == '+':
					if j == 0:
						aux_p.append(0)
					else:
						aux_p.append(size_sum)
					positives.append(location)
					names_p.append(name)
					chrms_p.append(chrx)
					if feature.type == "CDS":
						locus_p.append(locus_tag)
						if cogs_dict != None:
							cogs_p.append(cogs_dict.get(locus_tag))
				else:
					if j == 0:
						aux_n.append(0)
					else:
						aux_n.append(size_sum)
					negatives.append(location)
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
	#print(sizes_p)
	
	p_output = output + "_" + feat + "_pos.txt"
	n_output = output + "_" + feat + "_neg.txt"
	
	if divided:
		write_cog_files(new_pos, p_output, chrms_p, locus_p, cogs_p)
		write_cog_files(new_negs, n_output, chrms_n, locus_n, cogs_n)		
			 
	else:
		write_lines(new_pos, p_output, chrms_p, locus_p, cogs_p)
		write_lines(new_negs, n_output, chrms_n, locus_n, cogs_n)
	
	return


# Function to predict COG categories with DeepNOG
def get_categories(gbk_file, output, deepnog_confidence = 0):
	
	
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
	
	output_faa = output + ".faa"
	# command1 = "python " + path + "/genbank2faa.py --outputFile " + output_faa + " " + gbk_file
	
	try: 
		#print(command1)
		genbank2faa.genbankToFaa(gbk_file, output_faa)
		#process = subprocess.Popen(command1.split(), stdout=subprocess.PIPE)
		#output_, error = process.communicate()
		print()
		print("GBK file transformed into faa succesfully. File saved as", output_faa)
	except Exception as e:
		print(e)
		print("Error when transforming gbk to faa.")
		return
	
	# Predict from deepnog
	
	print("output", output)
	output_pred = output + "_prediction_deepnog.csv"
	command2 = "deepnog infer " + output_faa + " --out " + output_pred + " -db cog2020 -t 1"
	if deepnog_confidence > 0:
		command2 += " -c " + str(deepnog_confidence)
	
	try: 
		print()
		print("Deepnog prediction started")
		print()
		print(command2)
		process = subprocess.Popen(command2.split(), stdout=subprocess.PIPE)
		output_, error = process.communicate()
		print()
		print("Deepnog prediction finished succesfully. Predictions saved as", output_pred)
		print()
	except:
		print("Error when predicting categories from " + output_faa + " with deepnog.")
		return
	
	# Cross files 
	
	cogs_df = pd.read_csv(output_pred, sep=',', usecols=['sequence_id', 'prediction'])
	#tab_file = os.path.abspath(os.path.dirname(__file__)) + "/dataset/cog-20.def.tab" # To do: order this
	seek = []
	main = os.path.expanduser('~')
	dirs = list(map(str, listdir_r(main, "folder", seek)))
	if dirs == []:
		print("GenoVi folder not found")
		return
	dirs.sort(key=len)
	genovi_dir = dirs[0]

	tab_file =  genovi_dir + "/dataset/cog-20.def.tab"
	cogs_df.columns = ['id', 'cog']
	tab_df = pd.read_csv(tab_file, header=None, sep='\t', usecols=[0,1], encoding='cp1252') # (parche) We should check why this file is different if it's running on Windows. (Maybe it works on Linux too?)
	# tab_df = pd.read_csv(tab_file, header=None, sep='\t', usecols=[0,1]) # Original Linux form
	tab_df.columns = ['cog', 'category']
	merge_df = pd.merge(cogs_df, tab_df, on=["cog"])
	
	cogs_dict = merge_df.set_index('id')['category'].to_dict()
	
	return cogs_dict


# Base pipeline for complete genome.
def base_complete(gbk_file, output, cds, trna, get_cats, divided, k, init, end, sizes):
	
	flag = True
	
	if cds == False and get_cats == True:
		print("Error: Categories can only be predicted for CDS. Please enter output file paths for both CDS.") 
		flag = False
	elif get_cats == True:
		cogs_dict = get_categories(gbk_file, output)
	else:
		cogs_dict = None
		
	if divided == True and get_cats == False:
		print("Error: Division of categories is part of the categories prediction process. Include --get_categories in your command.") 
		flag = False
	
	if flag == True:
		
		create_kar_complete(output, k, init, end)
			
		if trna == True:
			create_feature_complete(gbk_file, output, sizes, k, "tRNA")
		if rrna == True:
			create_feature_complete(gbk_file, output, sizes, k, "rRNA")
		if cds == True:
			create_feature_complete(gbk_file, output, sizes, k, "CDS", cogs_dict)
		if divided == True:
			create_feature_complete(gbk_file, output, sizes, k, "CDS", cogs_dict, divided)

# Base pipeline for non-complete genome.	
def base(gbk_file, output, cds, trna, get_cats, divided, complete, rrna = False, deepnog_confidence = 0, verbose = False):
	
	flag = True

	cogs_p = set()
	cogs_n = set()
	
	if not cds and get_cats:
		print("Error: Categories can only be predicted for CDS. Please enter output file paths for both CDS.") 
		flag = False
	elif get_cats:
		cogs_dict = get_categories(gbk_file, output, deepnog_confidence)
	else:
		cogs_dict = None
		
	if divided and not get_cats:
		print("Error: Division of categories is part of the categories prediction process. Include --get_categories in your command.") 
		flag = False
	
	if flag:
		
		sizes, _, _ = create_kar(gbk_file, output, complete)
			
		if trna:
			create_feature(gbk_file, output, sizes, "tRNA", verbose = verbose)
		if rrna:
			create_feature(gbk_file, output, sizes, "rRNA", verbose = verbose)
		if cds:
			create_feature(gbk_file, output, sizes, "CDS", cogs_dict, verbose = verbose)
		if divided:
			cogs_p, cogs_n = create_feature(gbk_file, output, sizes, "CDS", cogs_dict, divided, verbose = verbose)
	return ((sizes, cogs_p, cogs_n))


if __name__ == '__main__':
	
	#gbk_file, output, cds_pos, cds_neg, trna_pos, trna_neg, get_cats, divided, complete = get_args()[:]
	gbk_file, output, cds, trna, rrna, get_cats, divided, complete = get_args()[:]
	
	try:
		gbk_name = gbk_file.split('/')[-1].split('.g')[-2]
	except:
		gbk_name = gbk_file.split('.g')[-2]
	
	if output == "":
		output = gbk_name + "/"
	elif output[-1] != "/":
		output = output + "/"
	
	if not os.path.isdir(output):
		os.mkdir(output)
	
	if complete == True:
		
		files_gbk = []
		gbk = open(gbk_file,"r")
		
		try:
			gbk_name = gbk_file.split("/")[-1].split(".")[-2]
		except:
			gbk_name = gbk_file.split(".")[-2]	
			
		sizes, inits, ends = create_kar(gbk_file, output, complete)

		for k, rec in enumerate(SeqIO.parse(gbk, "genbank")):
			
			output_ = output + 'replicon_' + str(k+1) + "/"
			
			if not os.path.isdir(output_):
				os.mkdir(output_)
			
			
			#filename = output_ + gbk_name + '_' + rec.id + '.gbk'
			filename = output_ + gbk_name + '_' + str(k+1) + '.gbk'
			files_gbk.append(filename)
			SeqIO.write([rec], open(filename, "w"), "genbank")
			
			output_ = output_ + gbk_name + '_' + str(k+1)
			
			base_complete(filename, output_, cds, trna, get_cats, divided, k, inits[k], ends[k], sizes)
	
	else:
		
		output_ = output + gbk_name
		base(gbk_file, output_, cds, trna, get_cats, divided, complete, rrna)

