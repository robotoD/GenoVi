#! /usr/bin/env python 
##########################################################################
#skewi.py calculates a single SkewI value for a given chromosome
#Copyright (C) 2020 Jennifer Lu, jlu26@jhmi.edu
#
#This file is part of SkewIT 
#
#SkewIT is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 3 of the license, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program; if not, see <http://www.gnu.org/licenses/>.

#####################################################################
#Jennifer Lu, jennifer.lu717@gmail.com
#02/17/2020
#
#This program calculates a Skew Index (SkewI) 
#for a given genome within the range 0 to 1. Higher SkewI values
#indicate a strong GC Skew signal while lower SkewI values
#indicate a potentially low GC Skew signal
#
#Given a multi-fasta file, the program will output one SkewI per sequence
######################################################################
import sys, os, argparse 
from time import gmtime
from time import strftime
from Bio import SeqIO
import numpy as np
#####################################################################
def usage(): 
    sys.stderr.write("\n #########################################################################################\n")
    sys.stderr.write(" ################################## USAGE: SKEWI.PY ######################################\n")
    sys.stderr.write(" ## > python skewi.py -i SEQ.FASTA -o SKEW.TXT\t\t\t\t\t\t##\n")
    sys.stderr.write(" ## \t -i SEQ.FASTA............fasta file (multi-fasta permitted)\t\t\t##\n")
    sys.stderr.write(" ## Optional Parameters:\t\t\t\t\t\t\t\t##\n")
    sys.stderr.write(" ## \t -o SKEW.TXT.............output file (see below)\t\t\t\t##\n")
    sys.stderr.write(" ## \t -k WINDOW_SIZE..........length of subsequences for which to calculate gc skew\t##\n")
    sys.stderr.write(" ## \t    .....................[default:20000]\t\t\t\t\t##\n" )
    sys.stderr.write(" ## \t -f FREQUENCY............number of bases between the start of each window\t##\n")
    sys.stderr.write(" ## \t    .....................[default: k == f]\t\t\t\t\t##\n")
    sys.stderr.write(" ## \t --min-seq-len LEN.......set a minimum sequence length\t\t\t\t##\n")
    sys.stderr.write(" ## \t    .....................[default: 500,000bp]\t\t\t\t\t##\n")
    sys.stderr.write(" ## \t --complete/--all........only analyze sequences with 'complete' in header\t##\n")
    sys.stderr.write(" ## \t    .....................[default: --complete]\t\t\t\t\t##\n")
    sys.stderr.write(" ## \t --plasmid/--no-plasmid..include/exclude plasmid sequences from analysis\t##\n")
    sys.stderr.write(" ## \t    .....................[default: --no-plasmid]\t\t\t\t##\n")
    sys.stderr.write(" ## Output file: If no output file is provided, SkewI will be printed to standard out\t##\n")
    sys.stderr.write(" ## \t      Otherwise, a tab-delimited file will be generated:\t\t\t##\n")
    sys.stderr.write(" ## \t      with two columns: 1) sequence ID 2) skewI value\t\t\t\t##\n")   
    sys.stderr.write(" #########################################################################################\n")
    sys.stderr.write(" #########################################################################################\n\n")
    exit(0)
#####################################################################
def main(): 
    parser = argparse.ArgumentParser()
    #Required parameter: 
    parser.add_argument("-i","--input","-s","--seq", 
        dest="in_file", required=False, default="",
        help="Sequence file for which to calculate gc_skew")
    #Defaults =20K/Same as window size
    parser.add_argument("-k","-w","--window-len",
        dest="window_size", required=False, default=20000, type=int,
        help="Window size for which to calculate each sign(g-c) [default: 20kb]")
    parser.add_argument("-f","--freq", dest="freq", type=int,  default=-1,
        help="Length between the start indices of each window. \
        [default: same as window_len, must be less than window size]")
    #No output provided = standard out 
    parser.add_argument('-o','--output', required=False, 
        dest="out_file",default="",
        help="Output text file to save skewi (skew index values)")
    #Complete Genomes/Plasmid Sequence options 
    parser.add_argument('--complete', required=False, dest='complete_tf',
        action='store_true', default=True, help='Only analyze complete genomes') 
    parser.add_argument('--all', required=False, dest='complete_tf',
        action='store_false', default=True, help='Analyze complete and draft genomes')
    parser.add_argument('--plasmid', required=False, dest='plasmid_tf',
        action='store_true', default=False, help='Analyze plasmid sequences') 
    parser.add_argument('--no-plasmid', required=False, dest='plasmid_tf',  
        action='store_false', default=False, help='Exclude plasmid sequences')
    #Sequence length minimum 
    parser.add_argument('--min-len', required=False, dest='min_seq_len',
        type=int, default=500000,
        help='Minimum sequence length to analyze [default: 500kb]')
    #Usage option
    parser.add_argument('--usage', required=False, dest='usage', 
        action='store_true', default=False,
        help='Prints usage information for this program')
    args=parser.parse_args() 
    
    ########################
    #Test parameters
    if args.usage:
        usage()
    if (args.in_file == ""):
        sys.stderr.write(" >> Please provide an input sequence file\n")
        usage()
        exit(1)
    if (not os.path.isfile(args.in_file)):
        sys.stderr.write(" >> ERROR: %s is not a valid file\n" % args.in_file)
        exit(1)
    #Lengths
    if (args.freq == -1):
        args.freq = args.window_size
    if (args.window_size < 0): 
        sys.stderr.write(" >> ERROR: --window-size must be greater than 0\n")
        exit(1)
    if (args.freq < 0):
        sys.stderr.write(" >> ERROR: --freq must be greater than 0\n")
        exit(1)
    if (args.freq > args.window_size):
        sys.stderr.write(" >> ERROR: window size must be >= frequency\n")
        exit(1) 
    ########################
    #Start program
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stderr.write(" >> PROGRAM START TIME: " + time + '\n')
    sys.stderr.write("    Input file: %s\n" % args.in_file) 
    if args.out_file != "":
        sys.stderr.write("    Output: %s\n" % args.out_file)
    else: 
        sys.stderr.write("    Output: system standard out\n")
    sys.stderr.write("    Window size (bp): %i\n" % args.window_size)
    sys.stderr.write("    Frequency (bp): %i\n" % args.freq) 
    sys.stderr.write("    Minimum sequence length (bp): %i\n" % args.min_seq_len)
    sys.stderr.write("    Complete sequences only: %s\n" % args.complete_tf)
    sys.stderr.write("    Exclude plasmids: %s\n" % (not args.plasmid_tf))
    sys.stderr.flush()
    ########################
    count_seqs = 0 
    argmax = []
    seq2skewi = {} 
    sys.stderr.write("\n")
    sys.stderr.write(" >> Processing file: %s\n" % args.in_file)
    sys.stderr.write("\r\t%i seqs evaluated: " % (count_seqs))
    for record in SeqIO.parse(args.in_file,'fasta'):
        my_description = str(record.description)
        if ("complete" not in my_description) and args.complete_tf:
            continue
        if ("plasmid" in my_description) and not args.plasmid_tf:
            continue
        my_seq = str(record.seq)
        if len(my_seq) < args.min_seq_len:
            continue
        count_seqs += 1
        #Get sequence and description
        #Print Update
        count = 0
        tot = len(my_seq)/args.window_size
        #Calculate skew
        skew = []
        for i in range(0,len(my_seq),args.window_size): 
            g = my_seq[i:i+args.window_size].count("G")
            c = my_seq[i:i+args.window_size].count("C")
            if (g-c) > 0:
                skew.append(1)
            elif (g-c) < 0:
                skew.append(-1)
            else:
                skew.append(0)
        #Final print 
        curr_size = len(skew)/2 
        sys.stdout.flush()
        skew += skew[:curr_size]
        #Calculate abs(T2-T1) where T is sum of values
        curr_diffs = []
        x = sum(skew[0:curr_size])
        y = sum(skew[curr_size:2*curr_size])
        for i in range(0, curr_size-1):
            curr_diffs.append(abs(x-y))
            x = x - skew[i] + skew[i+curr_size]
            y = y - skew[i+curr_size] + skew[i+2*curr_size]
        curr_diffs.append(abs(x-y))
        if len(curr_diffs) > 0:
            max_cd = float(max(curr_diffs))
            seq2skewi[my_description] = max_cd/float(len(my_seq))*float(args.window_size)
        sys.stderr.write("\r\t%i sequences evaluated" % (count_seqs))
        sys.stderr.flush()

    sys.stderr.write("\r\t%i total sequences evaluated (all finished)\n" % (count_seqs))
    sys.stderr.flush()

    if args.out_file == "":
        sys.stdout.write("Sequence\tSkewI\n")
        for seq in seq2skewi:
            sys.stdout.write("%s\t%0.10f\n" % (seq,seq2skewi[seq]))
    else:
        sys.stdout.write(" >> Printing SkewI values to file: %s\n" % args.out_file)
        o_file = open(args.out_file,'w')
        o_file.write("Sequence\tSkewI\n")
        for seq in seq2skewi:
            o_file.write("%s\t%0.10f\n" % (seq,seq2skewi[seq]))
        o_file.close()

    #End program
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stderr.write(" >> PROGRAM FINISH TIME: " + time + '\n')
if __name__== "__main__":
    main()
