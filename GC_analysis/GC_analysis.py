"""
A command-line utility for calculating the GC percentage and GC skew of a genomic sequence.
Based on Jennifer Lu's (jennifer.lu717@gmail.com) SkewIt and GC_analysis utility.
"""

import argparse as ap
import sys
from Bio import SeqIO

def get_args():
    parser = ap.ArgumentParser()
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-i", "--input_file", type=str, help="Name of the input file in FASTA format",
                               required=True)
    requiredNamed.add_argument("-w", "--window_size", type=int,
                               help="Number of base pairs where the GC percentage is calculated for",
                               default=1000)
    requiredNamed.add_argument("-s", "--shift", type=int, help="The shift increment. Defaults to window size", default=-1)
    parser.add_argument("-o", "--output_file", type=str, help="Name of the output files. Defaults to match input file.", default = " ")
    parser.add_argument("-ot", "--omit_tail", action="store_true", help="True: if the trailing sequence should be "
                                                                        "omitted. Default behaviour is to retain "
                                                                        "the leftover sequence.",
                        default=False)
    parser.add_argument("-f", "--output_format", type=str, choices=["wiggle"],
                        default="wiggle")
    parser.add_argument("-one", "--one_file", action="store_true", help="Force one file output", default=False)
    args = parser.parse_args()
    if args.shift == -1:
        args.shift = args.window_size
    if(args.output_file == " "):
        args.output_file = args.input_file
    return args.input_file, args.output_file, args.window_size, args.shift, args.omit_tail, args.output_format,\
           args.one_file

def write_content(loc, final_loc, data, file):
    file.write("chr" + str(globalIndex)
               + "\t" + str(loc + 1) + "\t" + str(final_loc) + "\t" + str(data) + "\n")

def generate_result():
    """
    Calculate GC percentage and write to output file.
    :return: None
    """
    seq_len = len(record)  # Get length of sequence
    i=0
    global max_GC_percentage
    global min_GC_percentage
    global maxSkew
    global minSkew
    i = -1
    for i in range((seq_len - window_size + shift) // shift):  # Iterate over the total number of shifts
        frag = record.seq[i * shift: i * shift + window_size]  # Extract the string for counting
        # Count number of C and G and convert to percentage
        percent = int(round((frag.count("C") + frag.count("G")) / float(window_size) * 100))
        g = float(frag.count("G"))
        c = float(frag.count("C"))
        if (g+c) > 0:
            new_calc = (g-c)/(g+c)
        else:
            new_calc = 0.0
        if(percent > max_GC_percentage): 
            max_GC_percentage = percent
        if(percent < min_GC_percentage):
            min_GC_percentage = percent
        if(new_calc > maxSkew): 
            maxSkew = new_calc
        if(new_calc < minSkew):
            minSkew = new_calc
        write_content(sequence_begin + i * shift, sequence_begin + i * shift + window_size, percent, result_content)
        write_content(sequence_begin + i * shift, sequence_begin + i * shift + window_size, new_calc, result_skew)
    if (i + 1) * shift < seq_len and not omit_tail:
        # if trailing sequence exits and omit_tail is False
        frag = record.seq[(i + 1) * shift:]
        percent = int(round((frag.count("C") + frag.count("G")) / float(len(frag)) * 100))
        g = float(frag.count("G"))
        c = float(frag.count("C"))
        if (g+c) > 0:
            new_calc = (g-c)/(g+c)
        else:
            new_calc = 0.0
        if(percent > max_GC_percentage): 
            max_GC_percentage = percent
        if(percent < min_GC_percentage):
            min_GC_percentage = percent
        if(new_calc > maxSkew): 
            maxSkew = new_calc
        if(new_calc < minSkew):
            minSkew = new_calc
        write_content(sequence_begin + (i + 1) * shift, sequence_begin + seq_len, percent, result_content)
        write_content(sequence_begin + (i + 1) * shift, sequence_begin + seq_len, new_calc, result_skew)
    return(seq_len)

if __name__ == "__main__":
    input_file, output_file, window_size, shift, omit_tail, output_format, one_file = get_args()[:]

    records = SeqIO.index(input_file, "fasta")
    records_num = len(records)
    sequence_begin = 0
    min_GC_percentage = 101
    max_GC_percentage = -1
    maxSkew = -float("inf")
    minSkew = float("inf")
    file = open(output_file + "_GC_content.wig", "w+")
    file2 = open(output_file + "_GC_skew.wig", "w+")
    file.close()
    file2.close()
    if records_num < 1:
        # No sequence in fasta file, corrupted
        sys.stdout.write("WARNING! {} contains no sequence data.\n".format(input_file))
        raise TypeError
    else:
        result_content = open(output_file + "_GC_content.wig", "a+")
        result_skew = open(output_file + "_GC_skew.wig", "a+")
        globalIndex = 1
        for record in SeqIO.parse(input_file, "fasta"):
            sequence_begin += generate_result()
            globalIndex += 1
        result_content.close()
        result_skew.close()
        file = open(output_file + ".maxmin", "w")
        file.write("#Plot GC content\nmin = {}\nmax = {}\n\n#Plot GC skew\nmin = {}\nmax = {}".format(min_GC_percentage, max_GC_percentage, minSkew, maxSkew))
        file.close()
