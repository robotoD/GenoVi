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
#
# This file defines an utility for transforming GenBank flat files into .faa FASTA files.
# Intended for GenoVi usage, might not work in other cases.

from Bio import SeqIO
import argparse
import re

# Modifies locus description in case that input file was generated
# with Prokka.
def modify_locus(input):
    input_handle  = open(input, "r")
    line = input_handle.readline()
    locus = ""
    locuses = []
    old = []
    while line:  
        if "LOCUS" in line:
            locus = line[:-1]
            old.append(locus)
            date = input_handle.readline()
            locus = locus + "   " + date
        elif "     source" in line:
            size = line.split("          ")[-1].split("..")[-1][:-1]
            if size in locus:
                aux = locus.split(size + " ")
                new_locus = aux[0] + "  " + size + " " + aux[1]
                locuses.append(new_locus)
        line = input_handle.readline()
    input_handle.close()
    
    input_modify  = open(input, "r")
    for i in range(len(locuses)):
        if i == 0:
            fileContent = re.sub(old[i], locuses[i], input_modify.read())
        else:
            fileContent = re.sub(old[i], locuses[i], fileContent)

    input_modify.close()
    input_modify = open(input, "w")
    input_modify.write(fileContent)
    input_modify.close()
    input_modify = open(input, "r")
        
# Function that transforms genbank file to faa file.   
def genbankToFaa(input, output, verbose = False):
    try:
        input_handle  = open(input, "r")
    except:
        input = input + ".gbk"
        input_handle  = open(input, "r")
    output_handle = open(output, "w")

    line = input_handle.readline()
    input_handle.close()
    if not line.split(' bp ')[-2].split(' ')[-1].isdecimal():
        modify_locus(input)
        
      
    for seq_record in SeqIO.parse(open(input,"r"), "genbank"):
        if verbose:
            print("Dealing with GenBank record %s" % seq_record.id)
        for seq_feature in seq_record.features :
            if seq_feature.type=="CDS":
                if "translation" in seq_feature.qualifiers and len(seq_feature.qualifiers['translation'])==1:
                    fasta = seq_feature.qualifiers['translation'][0]
                    output_handle.write(">%s from %s\n%s\n" % (
                        seq_feature.qualifiers['locus_tag'][0],
                        seq_record.name,
                        "\n".join([fasta[i:i+60] for i in range(0, len(fasta), 60)])))

    output_handle.close()
    input_handle.close()
    if verbose:
        print("Done")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Starting from a GenBank flat file, this simplifies it as a .fna")
    parser.add_argument("inputFile", help="flat file to be translated")
    parser.add_argument("--outputFile", help="output file name" )
    args = parser.parse_args()
    if ".gb" in args.inputFile:
        if args.inputFile[-4] == ".":
            args.inputFile = args.inputFile[:-4]
        elif args.inputFile[-3] == ".":
            args.inputFile = args.inputFile[:-3]
    if args.outputFile is None:
        args.outputFile = args.inputFile + "_converted.faa"
    genbankToFaa(args.inputFile, args.outputFile)
