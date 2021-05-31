#See my webpage:
#http://www.warwick.ac.uk/go/peter_cock/python/genbank2fasta/
from Bio import SeqIO
import argparse

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
    args.outputFile = args.inputFile + "_converted.fna"

try:
    input_handle  = open(args.inputFile, "r")
except:
    input_handle  = open(args.inputFile + ".gbk", "r")
output_handle = open(args.outputFile, "w")

#Short version:
#SeqIO.write(SeqIO.parse(input_handle, "genbank"), output_handle, "fasta")

#Long version, allows full control of fasta output
for seq_record in SeqIO.parse(input_handle, "genbank") :
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
print("Done")
