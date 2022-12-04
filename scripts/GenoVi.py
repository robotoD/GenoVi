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

#from create_raw import base, postprocess
#from GC_analysis import get_args_, write_content, generate_result, makeGC, createGC
#from genbank2fna import gbkToFna, mainFna 
#from genbank2faa import modify_locus, genbankToFaa, mainFaa
#from createConf import create_conf, create_conf_main
#from addText import addText 
#from mergeImages import mergeImages 
#from colours import parseColours
#from create_tables import gral_table, draw_histogram

from .addText import addText
from .colours import parseColours
from .create_raw import getArgs, listdir_r, ends_sorted, create_kar, create_feature, base_complete, base, get_categories, create_kar_complete, create_feature_complete, new_loc, write_cog_files, write_lines, createRaw, postprocess
from .createConf import create_conf, create_conf_main
from .GC_analysis import get_args_, write_content, generate_result, makeGC, createGC
from .genbank2faa import modify_locus, genbankToFaa, mainFaa
from .genbank2fna import gbkToFna, mainFna
from .mergeImages import mergeImages
from .create_tables import gral_table, draw_histogram
from scripts import __version__

import re
import argparse as ap
import os
from shutil import which
import pandas as pd
import gzip
import shutil


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
                outFile.write('<rect x="-1000" y="-1000" width="5000px" height="5000px" style="fill:{};"/>\n'.format(colour))
            else:
                outFile.write('<rect x="0" y="0" width="3000px" height="3000px" style="fill:{};"/>\n'.format(colour))
        else:
            outFile.write(line)
    if os.remove(fileName):
        raise RuntimeError("Theres been an error deleting a temporary file")
    if os.rename("Genovi_temp_file.svg", fileName):
        raise RuntimeError("Theres been an error deleting a temporary file")

def unzip(input_file):
    if input_file[-3:] == ".gz" or input_file[-2:] == ".z":
        old_input_file = input_file
        input_file = temp_folder + "/" + input_file[:-3 if input_file[-3] == "." else -2].split("/")[-1]
        with gzip.open(old_input_file, 'rb') as f_in:
            with open(input_file , 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)

# Full pipeline
# input: anotated genome filename.
def visualiseGenome(input_file, status, output_file = "genovi", 
                    cogs_classified = True, deepnog_confidence_threshold = 0, alignment = "center", scale = "variable", keep_temporary_files = False, reuse_predictions = False, window = 5000, verbose = False,
                    captions = True, captionsPosition = "auto", title = "", title_position = "center", italic_words = 2, size = False, tracks_explain = False,
                    colour_scheme = "auto", background_colour = "transparent", font_colour = "0, 0, 0", GC_content = "auto", GC_skew ='auto', tRNA = 'auto', rRNA = 'auto', CDS_positive = 'auto', CDS_negative = 'auto', skew_line_colour = '0, 0, 0',
                    wanted_cogs = "ABCDEFGHIJKLMNOPQRSTUVWXYZ+"):

    if not cairo:
        print("There's been an error finding cairoSVG library, so PNG images might be different from expected. Please prefer using SVG output.")
    if output_file[-4:] == ".svg" or output_file[-4:] == ".png":
        output_file = output_file[:-4]
    if deepnog_confidence_threshold > 1 or deepnog_confidence_threshold < 0:
        if verbose:
            print("DeepNOG lower bound must be between 0 and 1")
        raise Exception("DeepNOG lower bound must be between 0 and 1")
    wanted_cogs = wanted_cogs.upper()
    if wanted_cogs == "MET-" or wanted_cogs == "METABOLISM":
        wanted_cogs = ['C', 'E', 'F', 'G', 'H', 'I', 'P', 'Q']
    elif wanted_cogs == "CEL-"  or wanted_cogs == "CELLULAR PROCESSES AND SIGNALING":
        wanted_cogs = ['D', 'M', 'N', 'O', 'Y', 'U', 'V', 'W', 'Y', 'Z']
    elif wanted_cogs == "INF-" or wanted_cogs == "INFORMATION STORAGE AND PROCESSING":
        wanted_cogs = ['A', 'B', 'J', 'K', 'L', 'X']
    elif wanted_cogs == "POO-" or wanted_cogs == "POORLY CHARACTERIZED":
        wanted_cogs = ["R", "S", "None"]
    elif len(wanted_cogs) > 26:
        wanted_cogs = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z', 'None']
    else:
        try:
            wanted_cogs = int(wanted_cogs)
        except ValueError:
            wanted_cogs = list(wanted_cogs)
    

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
    if not os.path.exists(temp_folder) and not os.path.isdir(input_file):
        os.mkdir(temp_folder)
    if not os.path.exists(output_file):
        os.mkdir(output_file)
    unzip(input_file)
    if status == "complete":
        if os.path.isdir(input_file):
            input_folder = input_file
            os.chdir(output_file)
            output_folder = output_file
            files_to_draw = [x for x in os.listdir("../" + input_folder) if ".gb" in x]

            gral_stats_table = open("Gral_Stats.csv", "w")
            gral_stats_table.write("Replicon,Size (bp),GC content (%),CDS,tRNA,rRNA\n")
            if cogs_classified:
                cog_classification_table = open("COG_Classification.csv", "w")
                cog_classification_table.write(",Cellular Processes and Signaling,,,,,,,,,,Information Storage and Processing,,,,,,Metabolism,,,,,,,,Poorly Characterized,,,\nReplicon,D,M,N,O,T,U,V,W,Y,Z,A,B,J,K,L,X,C,E,F,G,H,I,P,Q,R,S,Unclassified\n")
            for input_file in files_to_draw:
                output_file = ".".join((input_file.split("/")[-1]).split(".")[0:-1])
                input_file = "../" + input_folder + "/" + input_file
                unzip(input_file)
                temp_folder = output_file + "-temp"
                if verbose:
                    print("processing file " + input_file)
                if not os.path.exists(temp_folder):
                    os.mkdir(temp_folder)
                if not os.path.exists(output_file):
                    os.mkdir(output_file)
                if input_file[-3:] == ".gz" or input_file[-2:] == ".z":
                    old_input_file = input_file
                    input_file = temp_folder + "/" + input_file[:-3 if input_file[-3] == "." else -2].split("/")[-1]
                    with gzip.open(old_input_file, 'rb') as f_in:
                        with open(input_file , 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)

                
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
                
                sizes_full = []
                lengths_full = []
                chrms_full = []
                gc_avg_full = []
                hists_full = []
                
                for i in range(1, len(contigs) + 1):
                    output_file_part = "contig_" + str(i) + "-" + output_file
                    file = temp_folder + "/" + str(i) + ".gbk"
                    if (not reuse_predictions) and os.path.exists(temp_folder + "/" + output_file_part + "_prediction_deepnog.csv"):
                        os.remove(temp_folder + "/contig_" + str(i) + "-" + output_file + "_prediction_deepnog.csv")
                    sizes, cogs_p, cogs_n, lengths, chrms, hist, wanted_cogs = base(file, temp_folder + "/" + output_file_part, output_file + "/" + output_file, True, True, cogs_classified, cogs_classified, False, True, deepnog_confidence_threshold, verbose, wanted_cogs=wanted_cogs)
                    #sizes_full = sizes_full + sizes
                    lengths_full = lengths_full + lengths
                    chrms_full = chrms_full + chrms

                    if hist is not None:
                        hist.columns = ["COG Category", "Frequency", "chr"+str(i)]
                        hists_full.append(hist)
                    
                    full_cogs = full_cogs.union(cogs_p).union(cogs_n)
                    images.append({"size": sizes[0], "fileName": output_file + "-contig_" + str(i) + ".svg"})
                    gbkToFna(file, temp_folder + "/" + output_file_part + ".fna", verbose)
                    maxmins, gc_avg = makeGC(temp_folder + "/" + output_file_part + ".fna", temp_folder + "/" + output_file_part, window)
                    gc_avg_full = gc_avg_full + gc_avg
                    
                    create_conf(output_file_part, temp_folder, maxmins, font_colour, GC_content, GC_skew, CDS_positive, CDS_negative, tRNA, rRNA, skew_line_colour, background_colour, cogs_classified, set(cogs_p).intersection(set(wanted_cogs)), set(cogs_n).intersection(set(wanted_cogs)), tracks_explain, 1)
                    
                    if verbose:
                        print("Drawing {}...".format(i))
                    os.system("circos circos.conf >/dev/null 2>&1")
                    os.system("circos -debug_group _all >/dev/null 2>&1")

                    if tracks_explain:
                        addText("", "center", "circos.svg", "circ_v.svg", captions = False, cogs_captions = True, font_colour = font_colour, tracks_explain = True)
                        os.remove("circos.svg")
                        os.rename("circ_v.svg", "circos.svg")

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
                
                chrms_full = postprocess(chrms_full)
                
                if len(hists_full) > 0:
                    new_hist = pd.DataFrame()
                    for hist in hists_full:
                        new_hist[hist.columns[-1]] = hist.iloc[:,-1]
                    new_hist['Frequency'] = new_hist.sum(axis=1)
                    new_hist['COG Category'] = hists_full[0]['COG Category']
                    new_hist = new_hist[['COG Category', 'Frequency']+[c for c in new_hist.columns if c[:3] == "chr"]]
                    
                    draw_histogram(new_hist, output_file + "/" + output_file, status) 
                
                gral_table(lengths_full, gc_avg_full, chrms_full, output_file + "/" + output_file + "_Gral_Stats.csv")

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
                        addText("", inFile = output_file + ".svg", captions = captions, cogs_captions = cogs_classified, captionsPosition = captionsPosition, cogs = full_cogs.intersection(set(wanted_cogs)),
                                pCDS_colour = CDS_positive, nCDS_colour = CDS_negative, tRNA_colour = tRNA, rRNA_colour = rRNA, GC_content_colour = GC_content, font_colour = font_colour)
                    else:
                        mergeImages(images, outFile = output_file + ".svg", align = alignment, scale = scale, background_colour = "none" if delete_background else background_colour)
                        addText(title, position = title_position, inFile = output_file + "/" + output_file + ".svg", italic = italic_words, captions = captions, cogs_captions = cogs_classified, captionsPosition = captionsPosition, cogs = full_cogs.intersection(set(wanted_cogs)),
                                pCDS_colour = CDS_positive, nCDS_colour = CDS_negative, tRNA_colour = tRNA, rRNA_colour = rRNA, GC_content_colour = GC_content, font_colour = font_colour)
                    os.remove(output_file + ".svg")
                    os.rename("titled_" + output_file + ".svg", output_file + "/" + output_file + ".svg")
                else:
                    mergeImages(images, outFile = output_file + "/" + output_file + ".svg", align = alignment, scale = scale, background_colour = "none" if delete_background else background_colour)
                for i in range(1, len(contigs) + 1):
                    os.rename(output_file + "-contig_" + str(i) + ".svg", output_file + "/" + output_file + "-contig_" + str(i) + ".svg")
                    os.rename(output_file + "-contig_" + str(i) + ".png", output_file + "/" + output_file + "-contig_" + str(i) + ".png")
                
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

                if cairo:
                    if verbose:
                        print("Converting to png...")
                    file = open(output_file + "/" + output_file + ".svg")
                    svg2png(bytestring = file.read(), write_to = output_file + "/" + output_file + ".png")
                    file.close()

                gral_stats_table.write(output_file + "\n")
                with open(output_file + "/" + output_file + "_Gral_Stats.csv") as local_table:
                    gral_stats_table.write("\n".join(local_table.read().split("\n")[1:]))
                if cogs_classified:
                    cog_classification_table.write(output_file + "\n")
                    with open(output_file + "/" + output_file + "_COG_Classification.csv") as local_table:
                        cog_classification_table.write("\n".join(local_table.read().split("\n")[2:]))
            
            gral_stats_table.close()
            if cogs_classified:
                cog_classification_table.close()
            try:
                os.remove("circos.png")
                # os.remove("circos.svg")
            except:
                pass
            os.chdir("..")
            output_file = output_folder
            input_file = input_folder

        else: # if it is only one input file
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
            
            sizes_full = []
            lengths_full = []
            chrms_full = []
            gc_avg_full = []
            hists_full = []
            
            for i in range(1, len(contigs) + 1):
                output_file_part = "contig_" + str(i) + "-" + output_file
                file = temp_folder + "/" + str(i) + ".gbk"
                if (not reuse_predictions) and os.path.exists(temp_folder + "/" + output_file_part + "_prediction_deepnog.csv"):
                    os.remove(temp_folder + "/contig_" + str(i) + "-" + output_file + "_prediction_deepnog.csv")
                sizes, cogs_p, cogs_n, lengths, chrms, hist, wanted_cogs = base(file, temp_folder + "/" + output_file_part, output_file + "/" + output_file, True, True, cogs_classified, cogs_classified, False, True, deepnog_confidence_threshold, verbose, wanted_cogs=wanted_cogs)
                #sizes_full = sizes_full + sizes
                lengths_full = lengths_full + lengths
                chrms_full = chrms_full + chrms

                if hist is not None:
                    hist.columns = ["COG Category", "Frequency", "chr"+str(i)]
                    hists_full.append(hist)
                
                full_cogs = full_cogs.union(cogs_p).union(cogs_n)
                images.append({"size": sizes[0], "fileName": output_file + "-contig_" + str(i) + ".svg"})
                gbkToFna(file, temp_folder + "/" + output_file_part + ".fna", verbose)
                maxmins, gc_avg = makeGC(temp_folder + "/" + output_file_part + ".fna", temp_folder + "/" + output_file_part, window)
                gc_avg_full = gc_avg_full + gc_avg
                
                create_conf(output_file_part, temp_folder, maxmins, font_colour, GC_content, GC_skew, CDS_positive, CDS_negative, tRNA, rRNA, skew_line_colour, background_colour, cogs_classified, set(cogs_p).intersection(set(wanted_cogs)), set(cogs_n).intersection(set(wanted_cogs)), tracks_explain, 1)
                
                if verbose:
                    print("Drawing {}...".format(i))
                os.system("circos circos.conf >/dev/null 2>&1")
                os.system("circos -debug_group _all >/dev/null 2>&1")

                if tracks_explain:
                    addText("", "center", "circos.svg", "circ_v.svg", captions = False, cogs_captions = True, font_colour = font_colour, tracks_explain = True)
                    os.remove("circos.svg")
                    os.rename("circ_v.svg", "circos.svg")

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
            
            chrms_full = postprocess(chrms_full)
            
            if len(hists_full) > 0:
                new_hist = pd.DataFrame()
                for hist in hists_full:
                    new_hist[hist.columns[-1]] = hist.iloc[:,-1]
                new_hist['Frequency'] = new_hist.sum(axis=1)
                new_hist['COG Category'] = hists_full[0]['COG Category']
                new_hist = new_hist[['COG Category', 'Frequency']+[c for c in new_hist.columns if c[:3] == "chr"]]
                
                draw_histogram(new_hist, output_file + "/" + output_file, status) 
            
            gral_table(lengths_full, gc_avg_full, chrms_full, output_file + "/" + output_file + "_Gral_Stats.csv")

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
                    addText("", inFile = output_file + ".svg", captions = captions, cogs_captions = cogs_classified, captionsPosition = captionsPosition, cogs = full_cogs.intersection(set(wanted_cogs)),
                            pCDS_colour = CDS_positive, nCDS_colour = CDS_negative, tRNA_colour = tRNA, rRNA_colour = rRNA, GC_content_colour = GC_content, font_colour = font_colour)
                else:
                    mergeImages(images, outFile = output_file + ".svg", align = alignment, scale = scale, background_colour = "none" if delete_background else background_colour)
                    addText(title, position = title_position, inFile = output_file + "/" + output_file + ".svg", italic = italic_words, captions = captions, cogs_captions = cogs_classified, captionsPosition = captionsPosition, cogs = full_cogs.intersection(set(wanted_cogs)),
                            pCDS_colour = CDS_positive, nCDS_colour = CDS_negative, tRNA_colour = tRNA, rRNA_colour = rRNA, GC_content_colour = GC_content, font_colour = font_colour)
                os.remove(output_file + ".svg")
                os.rename("titled_" + output_file + ".svg", output_file + "/" + output_file + ".svg")
            else:
                mergeImages(images, outFile = output_file + "/" + output_file + ".svg", align = alignment, scale = scale, background_colour = "none" if delete_background else background_colour)
            for i in range(1, len(contigs) + 1):
                os.rename(output_file + "-contig_" + str(i) + ".svg", output_file + "/" + output_file + "-contig_" + str(i) + ".svg")
                os.rename(output_file + "-contig_" + str(i) + ".png", output_file + "/" + output_file + "-contig_" + str(i) + ".png")
    else: # if status == "draft"
        if os.path.isdir(input_file):
            input_folder = input_file
            os.chdir(output_file)
            output_folder = output_file
            images = []
            files_to_draw = [x for x in os.listdir("../" + input_folder) if ".gb" in x]

            gral_stats_table = open("Gral_Stats.csv", "w")
            gral_stats_table.write("Replicon,Size (bp),GC content (%),CDS,tRNA,rRNA\n")
            if cogs_classified:
                cog_classification_table = open("COG_Classification.csv", "w")
                cog_classification_table.write(",Cellular Processes and Signaling,,,,,,,,,,Information Storage and Processing,,,,,,Metabolism,,,,,,,,Poorly Characterized,,,\nReplicon,D,M,N,O,T,U,V,W,Y,Z,A,B,J,K,L,X,C,E,F,G,H,I,P,Q,R,S,Unclassified\n")
            for input_file in files_to_draw:
                output_file = ".".join((input_file.split("/")[-1]).split(".")[0:-1])
                input_file = "../" + input_folder + "/" + input_file
                unzip(input_file)
                temp_folder = output_file + "-temp"
                if verbose:
                    print("processing file " + input_file)
                if not os.path.exists(temp_folder):
                    os.mkdir(temp_folder)
                if not os.path.exists(output_file):
                    os.mkdir(output_file)
                if input_file[-3:] == ".gz" or input_file[-2:] == ".z":
                    old_input_file = input_file
                    input_file = temp_folder + "/" + input_file[:-3 if input_file[-3] == "." else -2].split("/")[-1]
                    with gzip.open(old_input_file, 'rb') as f_in:
                        with open(input_file , 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)

                
                if (not reuse_predictions) and os.path.exists(temp_folder + "/" + output_file + "_prediction_deepnog.csv"):
                    os.remove(temp_folder + "/" + output_file + "_prediction_deepnog.csv")
                sizes, cogs_p, cogs_n, lengths, chrms, hist, wanted_cogs = base(input_file, temp_folder + "/" + output_file, output_file + "/" + output_file, True, True, cogs_classified, cogs_classified, False, True, deepnog_confidence_threshold, verbose, wanted_cogs=wanted_cogs) 
                images.append({"size": sum(sizes), "fileName": output_file + "/" + output_file + ".svg"})
                if hist is not None:
                    draw_histogram(hist, output_file + "/" + output_file, status)

                cogs_p = set(map(lambda x : "None" if x == None else x[0], cogs_p))
                cogs_n = set(map(lambda x : "None" if x == None else x[0], cogs_n))
                gbkToFna(input_file, temp_folder + "/" + output_file + ".fna", verbose)
                maxmins, gc_avg = makeGC(temp_folder + "/" + output_file + ".fna", temp_folder + "/" + output_file, window)
                create_conf(output_file, temp_folder, maxmins, font_colour, GC_content, GC_skew, CDS_positive, CDS_negative, tRNA, rRNA, skew_line_colour, background_colour, cogs_classified, cogs_p.intersection(set(wanted_cogs)), cogs_n.intersection(set(wanted_cogs)), tracks_explain, len(sizes))
                gral_table(lengths, gc_avg, chrms, output_file + "/" + output_file + "_Gral_Stats.csv")

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
                    if title == "filename":
                        addText(".".join(input_file.split("/")[-1].split(".")[:-1]), inFile = "circos.svg", captions = False, cogs_captions = False, size = sum(sizes) if size else "", tracks_explain = tracks_explain)
                    else:
                        addText("", inFile = "circos.svg", captions = False, cogs_captions = False, size = sum(sizes) if size else "", tracks_explain = tracks_explain)
                    os.remove("circos.svg")
                    os.rename("titled_circos.svg", output_file + "/" + output_file + ".svg")
                else:
                    os.rename("circos.svg", output_file + "/" + output_file + ".svg")
                    os.rename("circos.png", output_file + "/" + output_file + ".png")
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

                if cairo:
                    if verbose:
                        print("Converting to png...")
                    file = open(output_file + "/" + output_file + ".svg")
                    svg2png(bytestring = file.read(), write_to = output_file + "/" + output_file + ".png")
                    file.close()
                
                gral_stats_table.write(output_file + "\n")
                with open(output_file + "/" + output_file + "_Gral_Stats.csv") as local_table:
                    gral_stats_table.write("\n".join(local_table.read().split("\n")[1:]))
                if cogs_classified:
                    cog_classification_table.write(output_file + "\n")
                    with open(output_file + "/" + output_file + "_COG_Classification.csv") as local_table:
                        cog_classification_table.write("\n".join(local_table.read().split("\n")[2:]))
            
            gral_stats_table.close()
            if cogs_classified:
                cog_classification_table.close()
            mergeImages(images, outFile = output_folder + ".svg", align = alignment, scale = scale, background_colour = "none" if delete_background else background_colour)
            addText("", inFile=output_folder + ".svg", outFile="default", captions=captions, cogs_captions=cogs_classified, captionsPosition = captionsPosition,
                pCDS_colour = CDS_positive, nCDS_colour = CDS_negative, tRNA_colour = tRNA, rRNA_colour = rRNA, GC_content_colour = GC_content, font_colour = font_colour)
            os.remove(output_folder + ".svg")
            os.rename("titled_" + output_folder + ".svg", output_folder + ".svg")
            try:
                os.remove("circos.png")
            except:
                pass
            os.chdir("..")
            output_file = output_folder
            input_file = input_folder
        else:
            if (not reuse_predictions) and os.path.exists(temp_folder + "/" + output_file + "_prediction_deepnog.csv"):
                os.remove(temp_folder + "/" + output_file + "_prediction_deepnog.csv")
            sizes, cogs_p, cogs_n, lengths, chrms, hist, wanted_cogs = base(input_file, temp_folder + "/" + output_file, output_file + "/" + output_file, True, True, cogs_classified, cogs_classified, False, True, deepnog_confidence_threshold, verbose, wanted_cogs=wanted_cogs) 
                
            if hist is not None:
                draw_histogram(hist, output_file + "/" + output_file, status)

            cogs_p = set(map(lambda x : "None" if x == None else x[0], cogs_p))
            cogs_n = set(map(lambda x : "None" if x == None else x[0], cogs_n))
            gbkToFna(input_file, temp_folder + "/" + output_file + ".fna", verbose)
            maxmins, gc_avg = makeGC(temp_folder + "/" + output_file + ".fna", temp_folder + "/" + output_file, window)
            create_conf(output_file, temp_folder, maxmins, font_colour, GC_content, GC_skew, CDS_positive, CDS_negative, tRNA, rRNA, skew_line_colour, background_colour, cogs_classified, cogs_p.intersection(set(wanted_cogs)), cogs_n.intersection(set(wanted_cogs)), tracks_explain, len(sizes))
            gral_table(lengths, gc_avg, chrms, output_file + "/" + output_file + "_Gral_Stats.csv")

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
                addText(title, position = title_position, inFile = "circos.svg", italic = italic_words, captions = captions, cogs_captions = cogs_classified, captionsPosition = captionsPosition, cogs = cogs_p.union(cogs_n).intersection(set(wanted_cogs)),
                        pCDS_colour = CDS_positive, nCDS_colour = CDS_negative, tRNA_colour = tRNA, rRNA_colour = rRNA, GC_content_colour = GC_content, font_colour = font_colour, size = sum(sizes) if size else "", tracks_explain = tracks_explain)
                os.remove("circos.svg")
                os.rename("titled_circos.svg", output_file + "/" + output_file + ".svg")
            else:
                os.rename("circos.svg", output_file + "/" + output_file + ".svg")
                os.rename("circos.png", output_file + "/" + output_file + ".png")

    if cairo and not (status == "complete" and os.path.isdir(input_file)):
        if verbose:
            print("Converting to png...")
        file = open(output_file + "/" + output_file + ".svg")
        svg2png(bytestring = file.read(), write_to = output_file + "/" + output_file + ".png")
        file.close()
        if verbose and title != "" and italic_words > 0:
            print("\nWARNING! PNG version of image may look weird because italic text transformation is not yet properly implemented.\nPlease prefer using the svg version instead.\n")

    try:
        os.remove("circos.png")
    except:
        pass
    if (not keep_temporary_files) and (not os.path.isdir(input_file)):
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
    parser.add_argument("-i", "--input_file", type=str, help="Genbank file (.gb, .gbk, .gbff, .gbff.gz) path", required=True)
    parser.add_argument("-s", "--status", type=str, choices=["complete", "draft"], help="To draw each sequence as a unique circular representation (complete) or as a circle with bands per sequence (draft).", required = True)
    parser.add_argument("-o", "--output_file", type=str, help="Directory for output files. Default: genovi", default = "genovi")
    parser.add_argument("-cu", "--cogs_unclassified", action='store_false', help="Do not classify each protein sequence into COG categories.", required = False)
    parser.add_argument("-b", "--deepnog_confidence_threshold", type=float, help="Lower threshold for DeepNOG prediction certainty to be considered. Values in range [0,1] Default: 0", default = 0)
    parser.add_argument("--cogs", type=str, help="Symbol of each COG to draw. For example, for drawing only information storage and processing related, use 'ABJKLX'. When using a number n, the n most common categories will be drawn. By default, draws all of them.", default = "ABCDEFGHIJKLMNOPQRSTUVWXYZ+")
    parser.add_argument("-a", "--alignment", type=str, choices=["center", "top", "bottom", "A", "<", "U", "matrix", "two_lines"], help="When using --status complete, this defines the vertical alignment of every circular representation. Options: center, top, bottom, A (First on top), < (first to the left), U (Two on top, the rest below), matrix (multiple rows). By default this is defined by contig sizes", default = "auto")
    parser.add_argument("--scale", type=str, choices=["variable", "linear", "sqrt"], help="To select the scale-up ratio between each circular representations when the file is processes as a complete genome. This is useful to ensure visibility of each representation when the length difference is too high. Options: variable, linear, sqrt. Default: sqrt", default = "sqrt")
    parser.add_argument("-k", "--keep_temporary_files", action='store_true', help="Do not delete files used for circos image generation, including protein categories prediction by Deepnog.", required = False)
    parser.add_argument("-r", "--reuse_predictions", action='store_true', help="If available, reuse DeepNog prediction result from previous run. Useful only after --keep_temporary_files flag enabled.", required = False)
    
    parser.add_argument("-w", "--window", "--step", type=int, help="Base pair window for GC plotting. Default: 5000", default = 5000)
    parser.add_argument("-v", "--verbose", type=bool, help="Weather to print progress log.", default = True)

    text_group = parser.add_argument_group("text")
    text_group.add_argument("-c", "--captions_not_included", action='store_false', help = "Do not include captions for genomic features and COG colors.", required = False)
    text_group.add_argument("-cp", "--captions_position", type=str, choices=["right", "left", "top", "bottom", "auto"], help = "Where to position figure captions. Options: left, right, top, bottom, auto.", default = "auto")
    text_group.add_argument("-t", "--title", type=str, help="Title of the image (e.g., strain taxomomical identification). By default, it does not include title", default = "")
    text_group.add_argument("--title_position", type=str, choices=["center", "top", "bottom"], default = "center")
    text_group.add_argument("--italic_words", type=int, help="How many of words of the title should be written in italics. WARNING: Italics is not supported by the svg-to-png engine, so PNG image will not be available if this option is used. Default: 0", default = 0)
    text_group.add_argument("--size", action='store_true', help="To write the size (in base pairs) inside each circle.", required = False)
    text_group.add_argument("-te", "--tracks_explain", action='store_true', help="To include an additional text on each track to specify the meaning for each one.", required = False)

    colour_group = parser.add_argument_group("colours")
    colour_group.add_argument("-cs", "--colour_scheme", "--colour", type=str, help='''Colour scheme to use. Individual colours may be overriden wih other arguments. COG colours can not be edited.
                                Options: neutral, blue, purple, soil, grayscale, velvet, pastel, ocean, wood, beach, desert, ice, island, forest, toxic, fire, spring''', default = 'auto')
    colour_group.add_argument("-bc", "--background", "--background_colour", type=str, help="Background colour. Default: transparent", default = 'transparent')
    colour_group.add_argument("-fc", "--font_colour", type=str, help="Colour for ticks and caption texts. Default: black", default = '0, 0, 0')
    colour_group.add_argument("-pc", "--CDS_positive_colour", type=str, help="Colour for positive CDSs, in R, G, B format. Default: '180, 205, 222'", default = 'auto')
    colour_group.add_argument("-nc", "--CDS_negative_colour", type=str, help="Colour for negative CDSs, in R, G, B format. Default: '53, 176, 42'", default = 'auto')
    colour_group.add_argument("-tc", "--tRNA_colour", type=str, help="Colour for tRNAs, in R, G, B format. Default: '150, 5, 50'", default = 'auto')
    colour_group.add_argument("-rc", "--rRNA_colour", type=str, help="Colour for rRNAs, in R, G, B format. Default: '150, 150, 50'", default = 'auto')
    colour_group.add_argument("-cc", "--GC_content_colour", type=str, help="Colour for GC content, in R, G, B format. Default: '23, 0, 115'", default = "auto")
    colour_group.add_argument("-sc", "--GC_skew_colour", type=str, help="Colour scheme for GC skew. Might be a pair of RGB colours or Circos-readable code. For details, please read CIRCOS documentation. Default: '140, 150, 198 - 158, 188, 218'", default = 'auto')
    colour_group.add_argument("-sl", "--GC_skew_line_colour", type=str, help="Colour for GC skew line. Default: black", default = 'auto')
    
    parser.add_argument("--version", action="version", version=f'%(prog)s {__version__}')
    
    args = parser.parse_args()

    return (args.input_file, args.status, args.output_file,
    args.cogs_unclassified, args.deepnog_confidence_threshold, args.alignment, args.scale, args.keep_temporary_files, args.reuse_predictions, args.window, args.verbose,
    args.captions_not_included, args.captions_position, args.title, args.title_position, args.italic_words, args.size, args.tracks_explain,
    args.colour_scheme, args.background, args.font_colour, args.GC_content_colour, args.GC_skew_colour, args.tRNA_colour, args.rRNA_colour, args.CDS_positive_colour, args.CDS_negative_colour, args.GC_skew_line_colour,
    args.cogs)

def main():
    visualiseGenome(*get_args())


if __name__ == "__main__":
    visualiseGenome(*get_args())
