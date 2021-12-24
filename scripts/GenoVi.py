#!/usr/bin/python

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

#import create_raw as create_raw
#import GC_analysis as GC_analysis
#import genbank2fna as gbk2fna
#import genbank2faa as genbank2faa
#import createConf as createConf
#import addText as addText
#import mergeImages as merge
#import colors as colors

from .addText import addText
from .colors import parseColors
from .create_raw import getArgs, listdir_r, ends_sorted, create_kar, create_feature, base_complete, base, get_categories, create_kar_complete, create_feature_complete, new_loc, write_cog_files, write_lines, createRaw
from .createConf import create_conf, create_conf_main
from .GC_analysis import get_args_, write_content, generate_result, makeGC, createGC
from .genbank2faa import modify_locus, genbankToFaa, mainFaa
from .genbank2fna import gbkToFna, mainFna
from .mergeImages import mergeImages
from scripts import __version__

import re
import argparse as ap
import os
from shutil import which


__all__ = ['change_background', 'visualizeGenome', 'get_args', 'main',
           ]

try:
    from cairosvg import svg2png
    cairo = True
except:
    cairo = False
    # print("There's been an error finding cairoSVG library, so PNG images might be different from expected. Please prefer using SVG output.")

# Changes background color to a white Circos-generated svg image (modifies original file)
# input: Final color
def change_background(color, fileName = "circos.svg"):
    file = open(fileName)
    outFile = open("Genovi_temp_file.svg", "w")
    for line in file:
        if re.match('<rect x="-?1?0?0?0" y="0" width="[345]000px" height="3000px" style="fill:rgb\(255,255,255\);"/>',line):
            outFile.write('<rect x="-1000" y="-1000" width="5000px" height="5000px" style="fill:{};"/>'.format(color))
        else:
            outFile.write(line)
    if os.remove(fileName):
        raise RuntimeError("Theres been an error deleting a temporary file")
    if os.rename("Genovi_temp_file.svg", fileName):
        raise RuntimeError("Theres been an error deleting a temporary file")

# Full pipeline
# input: anotated genome filename.
def visualizeGenome(input_file, status, output_file = "circos", 
                    cogs_unclassified = True, deepnog_confidence_threshold = 0, alignment = "center", scale = "variable", keep_temporary_files = False, window = 5000, verbose = False,
                    captions = True, captionsPosition = "auto", title = "", title_position = "center", italic_words = 2, size = False,
                    color_scheme = "auto", background_color = "transparent", font_color = "0, 0, 0", GC_content = "auto", GC_skew ='auto', tRNA = 'auto', rRNA = 'auto', CDS_positive = 'auto', CDS_negative = 'auto', skew_line_color = '0, 0, 0'):

    if not cairo:
        print("There's been an error finding cairoSVG library, so PNG images might be different from expected. Please prefer using SVG output.")
    if output_file[-4:] == ".svg" or output_file[-4:] == ".png":
        output_file = output_file[:-4]
    if deepnog_confidence_threshold > 1 or deepnog_confidence_threshold < 0:
        if verbose:
            print("DeepNOG lower bound must be between 0 and 1")
        raise Exception("DeepNOG lower bound must be between 0 and 1")

    if which("circos") == None:
        if verbose:
            print("Circos is not installed. please install for using GenoVi.")
        raise Exception("Circos is not installed. please install for using GenoVi.")
    
    color_scheme, background_color, GC_content, GC_skew, tRNA, rRNA, CDS_positive, CDS_negative, skew_line_color = parseColors(color_scheme, background_color, GC_content, GC_skew, tRNA, rRNA, CDS_positive, CDS_negative, skew_line_color)
    delete_background = False
    if background_color == "transparent" or background_color == "none" or background_color == "auto":
        delete_background = True
        background_color = "white"
    
    if not os.path.exists("temp"):
        os.mkdir("temp")
    if status == "complete":
        file = open(input_file)
        contigs = file.read().split("\n//\n")
        file.close()
        if len(contigs[-1]) < 20:   # If file ends with "//"
            contigs = contigs[:-1]
        for i in range(1, len(contigs) + 1):
            contigFile = open("temp/" + str(i) + ".gbk", "w")
            contigFile.write(contigs[i - 1])
            contigFile.close()
        
        images = []
        full_cogs = set([])
        for i in range(1, len(contigs) + 1):
            file = "temp/" + str(i) + ".gbk"
            sizes, cogs_p, cogs_n = base(file, "temp/", True, True, cogs_unclassified, cogs_unclassified, False, True, deepnog_confidence_threshold, verbose)
            full_cogs = full_cogs.union(cogs_p).union(cogs_n)
            images.append({"size": sizes[0], "fileName": output_file + "-contig_" + str(i) + ".svg"})
            gbkToFna(file, "temp/gbk_converted.fna", verbose)
            maxmins = makeGC("temp/gbk_converted.fna", "temp/GC", window)
            create_conf(maxmins, font_color, GC_content, GC_skew, CDS_positive, CDS_negative, tRNA, rRNA, skew_line_color, background_color, cogs_unclassified, cogs_p, cogs_n)
            
            if verbose:
                print("Drawing {}...".format(i))
            os.system("circos circos.conf >/dev/null 2>&1")
            os.system("circos -debug_group _all >/dev/null 2>&1")
            if delete_background:
                change_background("none")
                if cairo:
                    if verbose:
                        print("Converting to png...")
                    svgFile = open("circos.svg")
                    svg2png(bytestring = svgFile.read(), write_to = output_file + ".png")
                    
                    svgFile.close()

            if size:
                addText("", "center", "circos.svg", "circ_v.svg", captions = False, cogs_captions = False, size = sizes[0], font_color = font_color)
                os.rename("circ_v.svg", output_file + "-contig_" + str(i) + ".svg")
            else:
                os.rename("circos.svg", output_file + "-contig_" + str(i) + ".svg")
            os.rename("circos.png", output_file + "-contig_" + str(i) + ".png")
            if cogs_unclassified:
                os.rename("temp/_prediction_deepnog.csv", "temp/" + str(i) + "_prediction_deepnog.csv")
            for file_i in os.listdir("temp/"):
                if ".gbk" not in file_i and ".csv" not in file_i and not "contig-" in file_i:
                    os.rename("temp/" + file_i, "temp/contig-" + str(i) + "-" + file_i)
            os.remove(file)

        if captions or title != "":
            if captionsPosition == "auto":
                captionsPosition = "top-right" if alignment == "bottom" else "bottom-right"
            if title_position == "center":
                addText(title, position = title_position, inFile = output_file + "-contig_" + "1.svg", italic = italic_words, captions = False, font_color = font_color)
                os.remove(output_file + "-contig_" + "1.svg")
                os.rename("titled_" + output_file + "-contig_" + "1.svg", output_file + "-contig_" + "1.svg")
                mergeImages(images, outFile = output_file + ".svg", align = alignment, scale = scale, background_color = "none" if delete_background else background_color)
                addText("", inFile = output_file + ".svg", captions = captions, cogs_captions = cogs_unclassified, captionsPosition = captionsPosition, cogs = full_cogs,
                                pCDS_color = CDS_positive, nCDS_color = CDS_negative, tRNA_color = tRNA, rRNA_color = rRNA, GC_content_color = GC_content, font_color = font_color)
            else:
                mergeImages(images, outFile = output_file + ".svg", align = alignment, scale = scale, background_color = "none" if delete_background else background_color)
                print(captionsPosition)
                addText(title, position = title_position, inFile = output_file + ".svg", italic = italic_words, captions = captions, cogs_captions = cogs_unclassified, captionsPosition = captionsPosition, cogs = full_cogs,
                                pCDS_color = CDS_positive, nCDS_color = CDS_negative, tRNA_color = tRNA, rRNA_color = rRNA, GC_content_color = GC_content, font_color = font_color)
            os.remove(output_file + ".svg")
            os.rename("titled_" + output_file + ".svg", output_file + ".svg")
        else:
            mergeImages(images, outFile = output_file + ".svg", align = alignment, scale = scale, background_color = "none" if delete_background else background_color)
        if not os.path.exists(output_file):
            os.mkdir(output_file)
        for i in range(1, len(contigs) + 1):
            os.rename(output_file + "-contig_" + str(i) + ".svg", output_file + "/" + output_file + "-contig_" + str(i) + ".svg")
            os.rename(output_file + "-contig_" + str(i) + ".png", output_file + "/" + output_file + "-contig_" + str(i) + ".png")
    else:
        sizes, cogs_p, cogs_n = base(input_file, "temp/", True, True, cogs_unclassified, cogs_unclassified, False, True, deepnog_confidence_threshold, verbose)
        cogs_p = set(map(lambda x : "None" if x == None else x[0], cogs_p))
        cogs_n = set(map(lambda x : "None" if x == None else x[0], cogs_n))
        
        gbkToFna(input_file, "temp/gbk_converted.fna", verbose)
        maxmins = makeGC("temp/gbk_converted.fna", "temp/GC", window)
        create_conf(maxmins, font_color, GC_content, GC_skew, CDS_positive, CDS_negative, tRNA, rRNA, skew_line_color, background_color, cogs_unclassified, cogs_p, cogs_n)

        if verbose:
            print("Drawing...")
        if which("circos") == None:
            if verbose:
                print("Circos is not installed. please install for using this program.")
            raise(Exception)
        os.system("circos circos.conf >/dev/null 2>&1")
        os.system("circos -debug_group _all >/dev/null 2>&1")
        if delete_background:
            change_background("none")
        if captions or title != "" or size:
            captionsPosition = "bottom-right" if captionsPosition == "auto" else captionsPosition
            addText(title, position = title_position, inFile = "circos.svg", italic = italic_words, captions = captions, cogs_captions = cogs_unclassified, captionsPosition = captionsPosition, cogs = cogs_p.union(cogs_n),
            pCDS_color = CDS_positive, nCDS_color = CDS_negative, tRNA_color = tRNA, rRNA_color = rRNA, GC_content_color = GC_content, font_color = font_color, size = sizes[0] if size else "")
            os.remove("circos.svg")
            os.rename("titled_circos.svg", "circos.svg")
    if cairo:
        if verbose:
            print("Converting to png...")
        file = open("circos.svg")
        svg2png(bytestring = file.read(), write_to = output_file + ".png")
        file.close()
    os.rename("circos.svg", output_file + ".svg")
    if output_file != "circos":
        os.remove("circos.png")

    if not keep_temporary_files:
        if verbose:
            print("deleting temporary files")
        os.remove("circos.conf")
        for file in os.listdir("temp/"):
            os.remove("temp/" + file)
        for file in os.listdir("conf/"):
            os.remove("conf/" + file)
        os.rmdir("temp")
        os.rmdir("conf")
    
def get_version():
	try:
		_dist = get_distribution('genovi')
		# Normalize case for Windows systems
		dist_loc = os.path.normcase(_dist.location)
		here = os.path.normcase(__file__)
		if not here.startswith(os.path.join(dist_loc, 'genovi')):
			# not installed, but there is another version that *is*
			raise DistributionNotFound
	except DistributionNotFound:
		__version__ = 'Please install this project with\n pip install genovi'
	else:
		__version__ = _dist.version
		
	return __version__
   
    
# Parse user arguments
def get_args():
    parser = ap.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str, help="Genbank file path", required=True)
    parser.add_argument("-s", "--status", type=str, choices=["complete", "draft"], help="To draw each contig as a complete circle by itself.", required = True)
    parser.add_argument("-o", "--output_file", type=str, help="Output image file path. Default: circos", default = "circos")
    parser.add_argument("-cu", "--cogs_unclassified", action='store_false', help="Do not clasify each coding sequence and draw them by color.", required = False)
    parser.add_argument("-b", "--deepnog_confidence_threshold", type=float, help="Lower bound for DeepNOG prediction certainty to be considered. Values in range [0,1] Default: 0", default = 0)
    parser.add_argument("-a", "--alignment", type=str, choices=["center", "top", "bottom", "A", "<", "U"], help="When using --status complete, this defines the vertical alignment of every contig. Options: center, top, bottom, A (First on top), < (first to the left), U (Two on top, the rest below). By default this is defined by contig sizes", default = "auto")
    parser.add_argument("--scale", type=str, choices=["variable", "linear", "sqrt"], help="When using --status complete, wether to use a different scale for tiny contigs, so to ensure visibility. Options: variable, linear, sqrt. Default: sqrt", default = "sqrt")
    parser.add_argument("-k", "--keep_temporary_files", action='store_true', help="Don't delete files used for circos image generation, including protein categories prediction by Deepnog.", required = False)
    parser.add_argument("-w", "--window", "--step", type=int, help="base pair window for CG plotting. Default: 5000", default = 5000)
    parser.add_argument("-v", "--verbose", type=bool, help="Wether to print progress logs.", default = True)

    text_group = parser.add_argument_group("text")
    text_group.add_argument("-c", "--captions_not_included", action='store_false', help = "Do not include color explanation.", required = False)
    text_group.add_argument("-cp", "--captions_position", type=str, choices=["right", "left", "auto"], help = "Where to insert color explanation. Options: left, right, auto.", default = "auto")
    text_group.add_argument("-t", "--title", type=str, help="Title of the image (strain name, or something like that). By default, it doesn't include title", default = "")
    text_group.add_argument("--title_position", type=str, choices=["center", "top", "bottom"], default = "center")
    text_group.add_argument("--italic_words", type=int, help="How many of the title's words should be written in italics. Default: 2", default = 2)
    text_group.add_argument("--size", action='store_true', help="To write the size (in base pairs) in each circle.", required = False)

    color_group = parser.add_argument_group("colors")
    color_group.add_argument("-cs", "--color_scheme", "--color", type=str, help='''Color scheme to use. Individual colors may be overriden wih other arguments. COGs' coloring can't be changed.
                                Options: neutral, blue, purple, soil, grayscale, velvet, pastel, ocean, wood, beach, desert, ice, island, forest, toxic, fire, spring''', default = 'auto')
    color_group.add_argument("-bc", "--background", "--background_color", type=str, help="Color for background. Default: transparent", default = 'transparent')
    color_group.add_argument("-fc", "--font_color", type=str, help="Color for ticks and caption texts. Default: black", default = '0, 0, 0')
    color_group.add_argument("-pc", "--CDS_positive_color", type=str, help="Color for positive CDSs, in R, G, B format. Default: '180, 205, 222'", default = 'auto')
    color_group.add_argument("-nc", "--CDS_negative_color", type=str, help="Color for negative CDSs, in R, G, B format. Default: '53, 176, 42'", default = 'auto')
    color_group.add_argument("-tc", "--tRNA_color", type=str, help="Color for tRNAs, in R, G, B format. Default: '150, 5, 50'", default = 'auto')
    color_group.add_argument("-rc", "--rRNA_color", type=str, help="Color for rRNAs, in R, G, B format. Default: '150, 150, 50'", default = 'auto')
    color_group.add_argument("-cc", "--GC_content_color", type=str, help="Color for GC content, in R, G, B format. Default: '23, 0, 115'", default = "auto")
    color_group.add_argument("-sc", "--GC_skew_color", type=str, help="Color scheme for GC skew. Might be a pair of RGB colors or Circos-understandable code. For details on this, please read CIRCOS documentation. Default: '140, 150, 198 - 158, 188, 218'", default = 'auto')
    color_group.add_argument("-sl", "--GC_skew_line_color", type=str, help="Color for GC skew line. Default: black", default = 'auto')
    
    parser.add_argument("--version", action="version", version=f'%(prog)s {__version__}')
    
    args = parser.parse_args()

    return (args.input_file, args.status, args.output_file,
    args.cogs_unclassified, args.deepnog_confidence_threshold, args.alignment, args.scale, args.keep_temporary_files, args.window, args.verbose,
    args.captions_not_included, args.captions_position, args.title, args.title_position, args.italic_words, args.size, 
    args.color_scheme, args.background, args.font_color, args.GC_content_color, args.GC_skew_color, args.tRNA_color, args.rRNA_color, args.CDS_positive_color, args.CDS_negative_color, args.GC_skew_line_color)

def main():
    visualizeGenome(*get_args())


if __name__ == "__main__":
    visualizeGenome(*get_args())
