"""
.. See the NOTICE file distributed with this work for additional information
   regarding copyright ownership.
   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
"""

"""
A command-line utility for calculating the GC percentage of a genomic sequence.
"""

import argparse as ap
import sys
from Bio import SeqIO


def get_args():
    """
    Helper function uses ArgParser to handler all the command-line input options.
    All arguments are keyword arguments. The required arguments are grouped under "required named arguments" section in
    help (-h).

    required named arguments:

    -i INPUT_FILE, --input_file INPUT_FILE
    INPUTFILE: Name of the input file in FASTA format

    -w WINDOW_SIZE, --window_size WINDOW_SIZE
    WINDOW_SIZE: Number of base pairs that the GC percentage is calculated for

    -s SHIFT, --shift SHIFT
    SHIFT: The shift increment (step size)

    optional arguments:

    -h, --help
    Show the help message and exit

    -o OUTPUT_FILE, --output_file OUTPUT_FILE
    OUTPUT_FILE: Name of the output file

    -ot, --omit_tail
    Use if the trailing sequence should be omitted. Default behaviour is to retain the leftover sequence.

    -f {wiggle}, --output_format {wiggle}
    Choose output formats from [wiggle] file.

    :returns: str, str/None, int, int, Bool, str
    """
    parser = ap.ArgumentParser()
    requiredNamed = parser.add_argument_group('required named arguments')
    requiredNamed.add_argument("-i", "--input_file", type=str, help="Name of the input file in FASTA format",
                               required=True)
    requiredNamed.add_argument("-w", "--window_size", type=int,
                               help="Number of base pairs where the GC percentage is calculated for",
                               required=True)
    requiredNamed.add_argument("-s", "--shift", type=int, help="The shift increment", default=-1)
    parser.add_argument("-o", "--output_file", type=str, help="Name of the output file")
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
    return args.input_file, args.output_file, args.window_size, args.shift, args.omit_tail, args.output_format,\
           args.one_file


def open_results_file():
    """
    A helper function to create the output file based on chosen output format and specified output filename.

    If output filename is given and output format is "wiggle", a text file called "OUTPUT_FILENAME.wig" will be opened
    with python's textIO wrapper, where OUTPUT_FILENAME is the string given after "-o" option in the commandline.

    If output filename is not specified, the result will be written to stdout following wiggle format.

    :return: file object
    """
    if output_file:
        file = open(output_file + ".wig", "w+")
    else:
        file = sys.stdout
    return file


def open_results_files():
    """
    A helper function to create the output file based on chosen output format and specified output filename. Used when
    multiple sequences present in one input file. "_seqNUM" will be added to the end of the output filename before the
    file extension.

    If output filename is given and output format is "wiggle", a text file called "OUTPUT_FILENAME_seqNUM.wig" will be
    opened with python's textIO wrapper, where OUTPUT_FILENAME is the string given after "-o" option in the commandline
    and NUM is the ordinal number of the sequence.

    If output filename is not specified, the result will be written to stdout following wiggle format.

    :return: file object
    """
    if output_file:
        file = open(output_file + ".wig", "a+")
    else:
        file = sys.stdout
    return file


def write_content(loc, final_loc, data):
    result.write("chr1\t" + str(loc + 1) + "\t" + str(final_loc) + "\t" + str(data) + "\n")


def generate_result():
    """
    Calculate GC percentage and write to output file.
    :return: None
    """
    seq_len = len(record)  # Get length of sequence
    i=0
    for i in range((seq_len - window_size + shift) // shift):  # Iterate over the total number of shifts
        frag = record.seq[i * shift: i * shift + window_size]  # Extract the string for counting
        # Count number of C and G and convert to percentage
        percent = int(round((frag.count("C") + frag.count("G")) / float(window_size) * 100))
        global max_GC_percentage
        if(percent > max_GC_percentage): 
            max_GC_percentage = percent
        global min_GC_percentage
        if(percent < min_GC_percentage):
            min_GC_percentage = percent
        write_content(sequence_begin + i * shift, sequence_begin + i * shift + window_size, percent)
    if (i + 1) * shift < seq_len and not omit_tail:
        # if trailing sequence exits and omit_tail is False
        frag = record.seq[(i + 1) * shift:]
        percent = int(round((frag.count("C") + frag.count("G")) / float(len(frag)) * 100))
        write_content(sequence_begin + (i + 1) * shift, sequence_begin + seq_len, percent)
    return(seq_len)


if __name__ == "__main__":
    error = []  # Store generated error message, and write to stderr at the end of stdout output
    input_file, output_file, window_size, shift, omit_tail, output_format, one_file = get_args()[:]
    new_output_format = output_format

    if output_format != "wiggle" and output_file is None:
        sys.stderr.write("WARNING! An output filename is needed to save output as {}. "
                         "The result is shown below:\n".format(output_format))
        error.append("WARNING! An output filename is needed to save output as {}. "
                     "The result is shown above.\n".format(output_format))
        new_output_format = "wiggle"

    output_format = new_output_format

    records = SeqIO.index(input_file, "fasta")
    records_num = len(records)
    sequence_begin = 0
    min_GC_percentage = 101
    max_GC_percentage = -1
    file = open(output_file + ".wig", "w+")
    file.close()
    if records_num < 1:
        # No sequence in fasta file, corrupted
        sys.stdout.write("WARNING! {} contains no sequence data.\n".format(input_file))
        raise TypeError
    elif records_num == 1 or one_file:
        result = open_results_file()
        # one sequence in fasta file or one output file for all sequences
        for record in SeqIO.parse(input_file, "fasta"):
            generate_result()
        result.close()
    else:
        result = open_results_files()
        for record in SeqIO.parse(input_file, "fasta"):
            sequence_begin += generate_result()
        result.close()
        file = open(output_file + ".maxmin", "w")
        file.write("max: {}\nmin: {}".format(max_GC_percentage, min_GC_percentage))
        file.close()
    if output_file is None:
        for err in error:
            sys.stderr.write(err)
