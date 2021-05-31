from Bio import GenBank
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
    args.outputFile = args.inputFile + "_converted.faa"

try:
    input_handle  = open(args.inputFile, "r")
except:
    input_handle  = open(args.inputFile + ".gbk", "r")
output_handle = open(args.outputFile, "w")

for seq_record in SeqIO.parse(input_handle, "genbank") :
    print("Dealing with GenBank record %s" % seq_record.id)
    for seq_feature in seq_record.features :
        if seq_feature.type=="CDS" :
            assert len(seq_feature.qualifiers['translation'])==1
            fasta = seq_feature.qualifiers['translation'][0]
            output_handle.write(">%s from %s\n%s\n" % (
                   seq_feature.qualifiers['locus_tag'][0],
                   seq_record.name,
                   "\n".join([fasta[i:i+60] for i in range(0, len(fasta), 60)])))

output_handle.close()
input_handle.close()
print("Done")
