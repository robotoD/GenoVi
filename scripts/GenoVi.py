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
#import colours as colours

from .addText import addText
from .colours import parseColours
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


__all__ = ['change_background', 'visualiseGenome', 'get_args', 'main',
           ]

try:
    from cairosvg import svg2png
    cairo = True
except:
    cairo = False
    # print("There's been an error finding cairoSVG library, so PNG images might be different from expected. Please prefer using SVG output.")

# Changes background colour to a white Circos-generated svg image (modifies original file)
# input: Final colour
def change_background(colour, finalImage = True, fileName = "circos.svg"):
    file = open(fileName)
    outFile = open("Genovi_temp_file.svg", "w")
    for line in file:
        if re.match('<rect (?:xmlns="http://www.w3.org/2000/svg" )?x="-?1?0?0?0" y="-?1?0?0?0" width="[345]000px" height="[345]000px" style="fill:rgb\([12]\d?\d,[12]\d?\d,[12]\d?\d\);"/>',line):
            if finalImage:
                outFile.write('<rect x="-1000" y="-1000" width="5000px" height="5000px" style="fill:{};"/>'.format(colour))
            else:
                outFile.write('<rect x="0" y="0" width="3000px" height="3000px" style="fill:{};"/>'.format(colour))
        else:
            outFile.write(line)
    if os.remove(fileName):
        raise RuntimeError("Theres been an error deleting a temporary file")
    if os.rename("Genovi_temp_file.svg", fileName):
        raise RuntimeError("Theres been an error deleting a temporary file")

# Full pipeline
# input: anotated genome filename.
def visualiseGenome(input_file, status, output_file = "circos", 
                    cogs_unclassified = True, deepnog_confidence_threshold = 0, alignment = "center", scale = "variable", keep_temporary_files = False, reuse_predictions = False, window = 5000, verbose = False,
                    captions = True, captionsPosition = "auto", title = "", title_position = "center", italic_words = 2, size = False,
                    colour_scheme = "auto", background_colour = "transparent", font_colour = "0, 0, 0", GC_content = "auto", GC_skew ='auto', tRNA = 'auto', rRNA = 'auto', CDS_positive = 'auto', CDS_negative = 'auto', skew_line_colour = '0, 0, 0'):

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
    
    colour_scheme, background_colour, GC_content, GC_skew, tRNA, rRNA, CDS_positive, CDS_negative, skew_line_colour = parseColours(colour_scheme, background_colour, GC_content, GC_skew, tRNA, rRNA, CDS_positive, CDS_negative, skew_line_colour)
    delete_background = False
    if background_colour == "transparent" or background_colour == "none" or background_colour == "auto":
        delete_background = True
        background_colour = "white"
    
    temp_folder = output_file + "-temp"
    if not os.path.exists(temp_folder):
        os.mkdir(temp_folder)
    if not os.path.exists(output_file):
        os.mkdir(output_file)
    if status == "complete":
        file = open(input_file)
        contigs = file.read().split("\n//\n")
        file.close()
        if len(contigs[-1]) < 20:   # If file ends with "//"
            contigs = contigs[:-1]
        for i in range(1, len(contigs) + 1):
            contigFile = open(temp_folder + "/" + str(i) + ".gbk", "w")
            contigFile.write(contigs[i - 1])
            contigFile.close()
        
        images = []
        full_cogs = set([])
        for i in range(1, len(contigs) + 1):
            output_file_part = "contig_" + str(i) + "-" + output_file
            file = temp_folder + "/" + str(i) + ".gbk"
            if (not reuse_predictions) and os.path.exists(temp_folder + "/" + output_file_part + "_prediction_deepnog.csv"):
                os.remove(temp_folder + "/contig_" + str(i) + "-" + output_file + "_prediction_deepnog.csv")
            sizes, cogs_p, cogs_n = base(file, temp_folder + "/" + output_file_part, True, True, cogs_unclassified, cogs_unclassified, False, True, deepnog_confidence_threshold, verbose)
            full_cogs = full_cogs.union(cogs_p).union(cogs_n)
            images.append({"size": sizes[0], "fileName": output_file + "-contig_" + str(i) + ".svg"})
            gbkToFna(file, temp_folder + "/" + output_file_part + ".fna", verbose)
            maxmins = makeGC(temp_folder + "/" + output_file_part + ".fna", temp_folder + "/" + output_file_part, window)
            create_conf(output_file_part, temp_folder, maxmins, font_colour, GC_content, GC_skew, CDS_positive, CDS_negative, tRNA, rRNA, skew_line_colour, background_colour, cogs_unclassified, cogs_p, cogs_n)
            
            if verbose:
                print("Drawing {}...".format(i))
            os.system("circos circos.conf >/dev/null 2>&1")
            os.system("circos -debug_group _all >/dev/null 2>&1")
            if delete_background:
                change_background("none", False)
                if cairo:
                    if verbose:
                        print("Converting to png...")
                    svgFile = open("circos.svg")
                    svg2png(bytestring = svgFile.read(), write_to = output_file + "-contig_" + str(i) + ".png")
                    svgFile.close()
                    os.remove("circos.png")

            if size:
                addText("", "center", "circos.svg", "circ_v.svg", captions = False, cogs_captions = False, size = sizes[0], font_colour = font_colour)
                os.rename("circ_v.svg", output_file + "-contig_" + str(i) + ".svg")
            else:
                os.rename("circos.svg", output_file + "-contig_" + str(i) + ".svg")
            if not delete_background:
                os.rename("circos.png", output_file + "-contig_" + str(i) + ".png")
            os.remove(file)

        if captions or title != "":
            if captionsPosition == "auto":
                if scale == "variable":
                    captionsPosition = "right"
                else:
                    captionsPosition = "top-right" if alignment == "bottom" else "bottom-right"
            if title_position == "center":
                addText(title, position = title_position, inFile = output_file + "-contig_" + "1.svg", italic = italic_words, captions = False, font_colour = font_colour)
                os.remove(output_file + "-contig_" + "1.svg")
                os.rename("titled_" + output_file + "-contig_" + "1.svg", output_file + "-contig_" + "1.svg")
                mergeImages(images, outFile = output_file + ".svg", align = alignment, scale = scale, background_colour = "none" if delete_background else background_colour)
                addText("", inFile = output_file + ".svg", captions = captions, cogs_captions = cogs_unclassified, captionsPosition = captionsPosition, cogs = full_cogs,
                                pCDS_colour = CDS_positive, nCDS_colour = CDS_negative, tRNA_colour = tRNA, rRNA_colour = rRNA, GC_content_colour = GC_content, font_colour = font_colour)
            else:
                mergeImages(images, outFile = output_file + ".svg", align = alignment, scale = scale, background_colour = "none" if delete_background else background_colour)
                print(captionsPosition)
                addText(title, position = title_position, inFile = output_file + "/" + output_file + ".svg", italic = italic_words, captions = captions, cogs_captions = cogs_unclassified, captionsPosition = captionsPosition, cogs = full_cogs,
                                pCDS_colour = CDS_positive, nCDS_colour = CDS_negative, tRNA_colour = tRNA, rRNA_colour = rRNA, GC_content_colour = GC_content, font_colour = font_colour)
            os.remove(output_file + ".svg")
            os.rename("titled_" + output_file + ".svg", output_file + "/" + output_file + ".svg")
        else:
            mergeImages(images, outFile = output_file + "/" + output_file + ".svg", align = alignment, scale = scale, background_colour = "none" if delete_background else background_colour)
        for i in range(1, len(contigs) + 1):
            os.rename(output_file + "-contig_" + str(i) + ".svg", output_file + "/" + output_file + "-contig_" + str(i) + ".svg")
            os.rename(output_file + "-contig_" + str(i) + ".png", output_file + "/" + output_file + "-contig_" + str(i) + ".png")
    else:
        if (not reuse_predictions) and os.path.exists(temp_folder + "/" + output_file + "_prediction_deepnog.csv"):
            os.remove(temp_folder + "/" + output_file + "_prediction_deepnog.csv")
        sizes, cogs_p, cogs_n = base(input_file, temp_folder + "/" + output_file, True, True, cogs_unclassified, cogs_unclassified, False, True, deepnog_confidence_threshold, verbose)
        cogs_p = set(map(lambda x : "None" if x == None else x[0], cogs_p))
        cogs_n = set(map(lambda x : "None" if x == None else x[0], cogs_n))
        
        gbkToFna(input_file, temp_folder + "/" + output_file + ".fna", verbose)
        maxmins = makeGC(temp_folder + "/" + output_file + ".fna", temp_folder + "/" + output_file, window)
        create_conf(output_file, temp_folder, maxmins, font_colour, GC_content, GC_skew, CDS_positive, CDS_negative, tRNA, rRNA, skew_line_colour, background_colour, cogs_unclassified, cogs_p, cogs_n)

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
            pCDS_colour = CDS_positive, nCDS_colour = CDS_negative, tRNA_colour = tRNA, rRNA_colour = rRNA, GC_content_colour = GC_content, font_colour = font_colour, size = sizes[0] if size else "")
            os.remove("circos.svg")
            os.rename("titled_circos.svg", output_file + "/" + output_file + ".svg")
        else:
            os.rename("circos.svg", output_file + "/" + output_file + ".svg")
            os.rename("circos.png", output_file + "/" + output_file + ".png")
    if cairo:
        if verbose:
            print("Converting to png...")
        file = open(output_file + "/" + output_file + ".svg")
        svg2png(bytestring = file.read(), write_to = output_file + "/" + output_file + ".png")
        file.close()

    try:
        os.remove("circos.png")
    except:
        pass

    if not keep_temporary_files:
        if verbose:
            print("deleting temporary files")
        os.remove("circos.conf")
        for file in os.listdir(temp_folder + "/"):
            if (not reuse_predictions) or "_prediction_deepnog.csv" != file[-23:]:
                os.remove(temp_folder + "/" + file)
        for file in os.listdir("conf/"):
            os.remove("conf/" + file)
        if not reuse_predictions:
            os.rmdir(temp_folder)
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
    parser.add_argument("-cu", "--cogs_unclassified", action='store_false', help="Do not clasify each coding sequence and draw them by colour.", required = False)
    parser.add_argument("-b", "--deepnog_confidence_threshold", type=float, help="Lower bound for DeepNOG prediction certainty to be considered. Values in range [0,1] Default: 0", default = 0)
    parser.add_argument("-a", "--alignment", type=str, choices=["center", "top", "bottom", "A", "<", "U"], help="When using --status complete, this defines the vertical alignment of every contig. Options: center, top, bottom, A (First on top), < (first to the left), U (Two on top, the rest below). By default this is defined by contig sizes", default = "auto")
    parser.add_argument("--scale", type=str, choices=["variable", "linear", "sqrt"], help="When using --status complete, wether to use a different scale for tiny contigs, so to ensure visibility. Options: variable, linear, sqrt. Default: sqrt", default = "sqrt")
    parser.add_argument("-k", "--keep_temporary_files", action='store_true', help="Don't delete files used for circos image generation, including protein categories prediction by Deepnog.", required = False)
    parser.add_argument("-r", "--reuse_predictions", action='store_true', help="If available, reuse DeepNog prediction result from previous run. This speeds up the process and is useful only after running this script with --keep_temporary_files flag enabled.", required = False)
    
    parser.add_argument("-w", "--window", "--step", type=int, help="base pair window for CG plotting. Default: 5000", default = 5000)
    parser.add_argument("-v", "--verbose", type=bool, help="Wether to print progress logs.", default = True)

    text_group = parser.add_argument_group("text")
    text_group.add_argument("-c", "--captions_not_included", action='store_false', help = "Do not include colour explanation.", required = False)
    text_group.add_argument("-cp", "--captions_position", type=str, choices=["right", "left", "top", "bottom", "auto"], help = "Where to insert colour explanation. Options: left, right, auto.", default = "auto")
    text_group.add_argument("-t", "--title", type=str, help="Title of the image (strain name, or something like that). By default, it doesn't include title", default = "")
    text_group.add_argument("--title_position", type=str, choices=["center", "top", "bottom"], default = "center")
    text_group.add_argument("--italic_words", type=int, help="How many of the title's words should be written in italics. Default: 2", default = 2)
    text_group.add_argument("--size", action='store_true', help="To write the size (in base pairs) in each circle.", required = False)

    colour_group = parser.add_argument_group("colours")
    colour_group.add_argument("-cs", "--colour_scheme", "--colour", type=str, help='''Colour scheme to use. Individual colours may be overriden wih other arguments. COGs' colouring can't be changed.
                                Options: neutral, blue, purple, soil, grayscale, velvet, pastel, ocean, wood, beach, desert, ice, island, forest, toxic, fire, spring''', default = 'auto')
    colour_group.add_argument("-bc", "--background", "--background_colour", type=str, help="Colour for background. Default: transparent", default = 'transparent')
    colour_group.add_argument("-fc", "--font_colour", type=str, help="Colour for ticks and caption texts. Default: black", default = '0, 0, 0')
    colour_group.add_argument("-pc", "--CDS_positive_colour", type=str, help="Colour for positive CDSs, in R, G, B format. Default: '180, 205, 222'", default = 'auto')
    colour_group.add_argument("-nc", "--CDS_negative_colour", type=str, help="Colour for negative CDSs, in R, G, B format. Default: '53, 176, 42'", default = 'auto')
    colour_group.add_argument("-tc", "--tRNA_colour", type=str, help="Colour for tRNAs, in R, G, B format. Default: '150, 5, 50'", default = 'auto')
    colour_group.add_argument("-rc", "--rRNA_colour", type=str, help="Colour for rRNAs, in R, G, B format. Default: '150, 150, 50'", default = 'auto')
    colour_group.add_argument("-cc", "--GC_content_colour", type=str, help="Colour for GC content, in R, G, B format. Default: '23, 0, 115'", default = "auto")
    colour_group.add_argument("-sc", "--GC_skew_colour", type=str, help="Colour scheme for GC skew. Might be a pair of RGB colours or Circos-understandable code. For details on this, please read CIRCOS documentation. Default: '140, 150, 198 - 158, 188, 218'", default = 'auto')
    colour_group.add_argument("-sl", "--GC_skew_line_colour", type=str, help="Colour for GC skew line. Default: black", default = 'auto')
    
    parser.add_argument("--version", action="version", version=f'%(prog)s {__version__}')
    
    args = parser.parse_args()

    return (args.input_file, args.status, args.output_file,
    args.cogs_unclassified, args.deepnog_confidence_threshold, args.alignment, args.scale, args.keep_temporary_files, args.reuse_predictions, args.window, args.verbose,
    args.captions_not_included, args.captions_position, args.title, args.title_position, args.italic_words, args.size, 
    args.colour_scheme, args.background, args.font_colour, args.GC_content_colour, args.GC_skew_colour, args.tRNA_colour, args.rRNA_colour, args.CDS_positive_colour, args.CDS_negative_colour, args.GC_skew_line_colour)

def main():
    visualiseGenome(*get_args())


if __name__ == "__main__":
    visualiseGenome(*get_args())
