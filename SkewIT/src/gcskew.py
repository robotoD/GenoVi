#! /usr/bin/env python 
##########################################################################
#gcskew.py calculates gcskew values for a given chromosome
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
##############################################################################
#Jennifer Lu, jennifer.lu717@gmail.com
#02/26/2020
#
#This program calculate gc_skew for a genome provided to this program.
#By default, the program will calculate gc-skew for 20kb windows every 20kb 
#(adjacent/non-overlapping windows of 20kb)
#
#Users can specify window size and frequency of the calculation. If a specified
#frequency is less than the window size, windows will be overlapping 
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
    sys.stderr.write(" ## > python gcskew.py -i SEQ.FASTA -o GCSKEW.TXT\t\t\t\t\t##\n")
    sys.stderr.write(" ## \t -i SEQ.FASTA............fasta file (multi-fasta permitted)\t\t\t##\n")
    sys.stderr.write(" ## Optional Parameters:\t\t\t\t\t\t\t\t##\n")
    sys.stderr.write(" ## \t -o SKEW.TXT.............output file (see below)\t\t\t\t##\n")
    sys.stderr.write(" ## \t -k WINDOW_SIZE..........length of subsequences for which to calculate gc skew\t##\n")
    sys.stderr.write(" ## \t    .....................[default:20000]\t\t\t\t\t##\n" )
    sys.stderr.write(" ## \t -f FREQUENCY............number of bases between the start of each window\t##\n")
    sys.stderr.write(" ## \t    .....................[default: k == f]\t\t\t\t\t##\n")
    sys.stderr.write(" ## Output file: If no output file is provided, SkewI will be printed to gcskew.txt\t##\n")
    sys.stderr.write(" ## \t      The output file is a 3 column, tab-delimited file listing:\t\t##\n")
    sys.stderr.write(" ## \t      (1) sequence ID (2) starting index of window (3) gc skew value\t\t##\n")   
    sys.stderr.write(" #########################################################################################\n")
    sys.stderr.write(" #########################################################################################\n\n")
    exit(0)
#####################################################################
def main(): 
    parser = argparse.ArgumentParser()
    parser.add_argument("-i","--input","-s","--seq", 
        dest="in_file", required=False, default="",
        help="FASTA/Multi-FASTA sequence file for which to calculate gc_skew")
    parser.add_argument("-k","-w","--window-size", default=20000,
        dest="window_size", required=False, type=int,
        help="Window size for which to calculate each g-c/g+c [default: 20Kbp]")
    parser.add_argument("-f","--freq",required=False,type=int,
        dest="freq", default=-1, 
        help="Frequency at which to calculate GC skew [default: same as window size]")
    parser.add_argument("-o","--output", required=False, default="gcskew.txt",
        dest="out_file", help="Output text file to save gc skew calculations [default: gcskew.txt]")
    parser.add_argument('--usage', required=False, dest='usage', 
        action='store_true', default=False,
        help='Prints usage information for this program')
    args=parser.parse_args() 
    
    ##############################################
    #Test Parameters
    if args.usage:
        usage()
    if (args.in_file == ""):
        sys.stderr.write(" >> Please provide an input sequence file\n")
        usage()
        exit(1)
    if (not os.path.isfile(args.in_file)):
        sys.stderr.write(" >> ERROR: %s is not a valid file\n" % args.in_file) 
        exit(1)
    if args.freq == -1:
        freq = args.window_size
    else:
        freq = args.freq
    if (args.window_size < 0):
        sys.stderr.write(" >> ERROR: window size -w must be greater than 0\n")
        exit(1)
    if (freq < 0):
        sys.stderr.write(" >> ERROR: -f/--freq must be greater than 0\n")
        exit(1)
    if (args.freq > args.window_size):
        sys.stderr.write(" >> ERROR: window size must be >= frequency\n")
        exit(1)
    ##############################################
    #Start program
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("PROGRAM START TIME: " + time + '\n')
    sys.stdout.write("    Input file: %s\n" % args.in_file) 
    sys.stdout.write("    Output: %s\n" % args.out_file)
    sys.stdout.write("    Window size (bp): %i\n" % args.window_size)
    sys.stdout.write("    Frequency (bp): %i\n" % freq) 
    sys.stdout.flush()
    ############################################################################
    #Process sequence file 
    id2string = {}
    id2name = {}
    count_seqs = 0 
    sys.stdout.write(">> Reading sequence file %s\n" % args.in_file)
    sys.stdout.write("\t%i seqs found" % (count_seqs))
    sys.stdout.flush()
    for record in SeqIO.parse(args.in_file,'fasta'):
        count_seqs += 1
        sys.stdout.write("\r\t%i seqs found" % (count_seqs))
        sys.stdout.flush()
        #Save string
        id2name[count_seqs] = record.id
        id2string[count_seqs] = str(record.seq)
    sys.stdout.write("\r\t%i seqs found\n" % (count_seqs))
    sys.stdout.flush()

    #Calculate and plot gc skew
    tot = count_seqs
    count_seqs = 0
    sys.stdout.write(">> Calculating/Saving GC Skew for %i seqs\n" % tot)
    sys.stdout.write("\tworking on sequence %i/%i" % (count_seqs,tot))
    sys.stdout.flush()
    #Prep output file
    o_file = open(args.out_file,'w')
    # o_file.write("Sequence\tIndex\tGC Skew (%ikb)\n" % (args.window_size/1000))
    #For each sequence, calculate gc skew and print
    thisSequenceBeginsAt = 0
    for i in range(1,tot+1): 
        #Update 
        count_seqs += 1
        sys.stdout.write("\r\tworking on sequence %i/%i" % (count_seqs,tot))
        sys.stdout.flush()
        #Calculate gc skew
        my_seq = id2string[i].upper()
        my_description = id2name[i]
        count = 0
        #Calculate skew
        g = 0.0
        c = 0.0
        curr_seq = ""
        for i in range(0, len(my_seq), freq):
            curr_seq = my_seq[i:min([i+args.window_size, len(my_seq)])]
            g = float(curr_seq.count("G"))
            c = float(curr_seq.count("C"))
            if (g+c) > 0:
                new_calc = (g-c)/(g+c)
            else:
                new_calc = 0.0
            #Print to file
            endOfThisSegment = thisSequenceBeginsAt + min([i + freq - 1, len(my_seq)])
            o_file.write("%s\t%i\t%i\t%0.8f\n" % (my_description, i + thisSequenceBeginsAt, endOfThisSegment, new_calc))
        thisSequenceBeginsAt += len(my_seq) + 1
    sys.stdout.write("\r\tFinished calculating gc skew for %i sequences\n" % tot)
    sys.stdout.flush()
    o_file.close()
    #End program
    time = strftime("%m-%d-%Y %H:%M:%S", gmtime())
    sys.stdout.write("PROGRAM FINISH TIME: " + time + '\n')

if __name__== "__main__":
    main()
