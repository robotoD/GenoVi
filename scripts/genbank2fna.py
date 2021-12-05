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
# This file defines an utility for transforming GenBank flat files into .fna FASTA files.
# Intended for GenoVi usage, might not work in other cases.
# Based on Peter Cock pipeline: http://www.warwick.ac.uk/go/peter_cock/python/genbank2fasta/


from Bio import SeqIO
import argparse


def gbkToFna(input, output = None, verbose = False):
    if output is None:
        output = input + "_converted.fna"
    try:
        input_handle  = open(input, "r")
    except:
        input_handle  = open(input + ".gbk", "r")
    output_handle = open(output, "w")

    #Short version:
    #SeqIO.write(SeqIO.parse(input_handle, "genbank"), output_handle, "fasta")

    #Long version, allows full control of fasta output
    for seq_record in SeqIO.parse(input_handle, "genbank"):
        if verbose:
            print("Dealing with GenBank record %s" % seq_record.id)
        try: # The old way, removed in Biopython 1.73
            fasta = seq_record.seq.tostring()
        except AttributeError: # The new way, needs Biopython 1.45 or later.
            fasta = str(seq_record.seq)
        
        output_handle.write(">%s %s\n%s\n" % (
            seq_record.id,
            seq_record.description,
            "\n".join([fasta[i:i+60] for i in range(0, len(fasta), 60)])))

    output_handle.close()
    input_handle.close()
    if verbose:
        print("Done")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Starting from a GenBank flat file, this simplifies it as a .fna")
    parser.add_argument("-i", "--inputFile", help = "flat file to be translated", required = True)
    parser.add_argument("-o", "--outputFile", help = "output file name" )
    args = parser.parse_args()
    if ".gb" in args.inputFile:
        if args.inputFile[-4] == ".":
            args.inputFile = args.inputFile[:-4]
        elif args.inputFile[-3] == ".":
            args.inputFile = args.inputFile[:-3]
    gbkToFna(args.inputFile, args.outputFile)
