import argparse as ap
import scripts.create_raw as create_raw
import scripts.GC_analysis as GC_analysis
import os
import scripts.genbank2fna as gbk2fna
import scripts.createConf as createConf
import scripts.addText as addText
import scripts.mergeImages as merge
from shutil import which
import re
#from cairosvg import svg2png # The import is actually lower in the script and it is not always strictly necessary, so this should stay commented.

def change_background(color):
    file = open("circos.svg")
    outFile = open("Genovi_temp_file.svg", "w")
    for line in file:
        if '<rect x="0" y="0" width="3000px" height="3000px" style="fill:rgb(255,255,255);"/>' in line:
            outFile.write('<rect x="0" y="0" width="3000px" height="3000px" style="fill:{};"/>'.format(color))
        else:
            outFile.write(line)
    os.remove("circos.svg")
    os.rename("Genovi_temp_file.svg", "circos.svg")


def visualizeGenome(gbk_file, output = "circos", 
                    cogs = True, deepnog_bound = 0, legend = True, separate = False, circles_alignment = "center", scale = "variable", keep_temp_files = False, window = 5000,
                    title = "", titlePos = "center", italic = 2,
                    color_scheme = "auto", background_color = "none", GC_content = "auto", GC_skew ='auto', tRNA = 'auto', rRNA = 'auto', CDS_positive = 'auto', CDS_negative = 'auto', skew_line_color = '0, 0, 0'):

    if re.match("^\s*[012]?\d?\d\s*,\s*[012]?\d?\d\s*,\s*[012]?\d?\d\s*$", background_color):
        background_color = "rgb(" + background_color + ")"

    color_scheme = color_scheme.lower()
    if color_scheme == "auto" or color_scheme == "neutral":
        GC_content = "94, 120, 145" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("bupu-7-seq-%d",remap_int(var(value),0,0,4,3)))' if GC_skew == "auto" else GC_skew
        tRNA = "67, 14, 110" if tRNA == "auto" else tRNA
        rRNA = "67, 110, 14" if rRNA == "auto" else rRNA
        CDS_positive = "186, 186, 186" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "140, 140, 140" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "171, 171, 171" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "blue" or color_scheme == "blues":
        GC_content = "14, 29, 130" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("rdbu-7-div-%d",remap_int(var(value),0,0,7,5)))' if GC_skew == "auto" else GC_skew
        tRNA = "99, 103, 186" if tRNA == "auto" else tRNA
        rRNA = "89, 123, 186" if rRNA == "auto" else rRNA
        CDS_positive = "191, 204, 217" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "171, 178, 217" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "163, 191, 217" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "purple" or color_scheme == "purples":
        GC_content = "57, 1, 120" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("puor-7-div-%d",remap_int(var(value),0,0,7,5)))' if GC_skew == "auto" else GC_skew
        tRNA = "108, 95, 156" if tRNA == "auto" else tRNA
        rRNA = "128, 85, 156" if rRNA == "auto" else rRNA
        CDS_positive = "196, 187, 237" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "151, 143, 191" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "163, 163, 217" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "soil":
        GC_content = "89, 60, 6" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("brbg-7-div-%d",remap_int(var(value),0,0,7,5)))' if GC_skew == "auto" else GC_skew
        tRNA = "0, 112, 11" if tRNA == "auto" else tRNA
        rRNA = "20, 112, 1" if rRNA == "auto" else rRNA
        CDS_positive = "185, 199, 186" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "175, 201, 178" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "89, 194, 115" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "greyscale" or color_scheme == "grayscale" or color_scheme == "grey" or color_scheme == "gray":
        GC_content = "87, 87, 87" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("rdgy-7-div-%d",remap_int(var(value),0,0,7,5)))' if GC_skew == "auto" else GC_skew
        tRNA = "115, 115, 115" if tRNA == "auto" else tRNA
        rRNA = "115, 115, 115" if rRNA == "auto" else rRNA
        CDS_positive = "209, 209, 209" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "184, 184, 184" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "171, 171, 171" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "velvet" or color_scheme == "pink" or color_scheme == "red":
        GC_content = "82, 1, 30" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("purd-7-seq-%d",remap_int(var(value),0,0,7,4)))' if GC_skew == "auto" else GC_skew
        tRNA = "130, 1, 49" if tRNA == "auto" else tRNA
        rRNA = "130, 21, 69" if rRNA == "auto" else rRNA
        CDS_positive = "217, 195, 205" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "196, 187, 193" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "230, 200, 211" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "pastel":
        GC_content = "196, 227, 255" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("pastel1-7-qual-%d",remap_int(var(value),0,0,4,3)))' if GC_skew == "auto" else GC_skew
        tRNA = "143, 143, 143" if tRNA == "auto" else tRNA
        rRNA = "143, 143, 143" if rRNA == "auto" else rRNA
        CDS_positive = "255, 222, 235" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "181, 255, 214" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "150, 150, 150" if skew_line_color == "auto" else skew_line_color
    
    if not os.path.exists("temp"):
        os.mkdir("temp")
    if separate:
        file = open(gbk_file)
        contigs = file.read().split("\n//\n")
        file.close()
        if len(contigs[-1]) < 20:   # Si el archivo termina con un "//"
            contigs = contigs[:-1]
        for i in range(1, len(contigs) + 1):
            contigFile = open("temp/" + str(i) + ".gbk", "w")
            contigFile.write(contigs[i - 1])
            contigFile.close()
        
        images = []
        i = 1
        full_cogs = set([])
        for contig in contigs:
            file = "temp/" + str(i) + ".gbk"
            sizes, cogs_p, cogs_n = create_raw.base(file, "temp/", True, True, cogs, cogs, False, True, deepnog_bound)
            full_cogs = full_cogs.union(cogs_p).union(cogs_n)
            images.append({"size": sizes[0], "fileName": str(i) + ".svg"})
            gbk2fna.gbkToFna(file, "temp/gbk_converted.fna")
            maxmins = GC_analysis.makeGC("temp/gbk_converted.fna", "temp/GC", window)
            createConf.create_conf(maxmins, GC_content, GC_skew, CDS_positive, CDS_negative, tRNA, rRNA, skew_line_color, cogs, cogs_p, cogs_n)
            
            print("Drawing {}...".format(i))
            if which("circos") == None:
                print("Circos is not installed. please install for using this program.")
                raise(Exception)
            os.system("circos circos.conf >/dev/null 2>&1")
            os.system("circos -debug_group _all >/dev/null 2>&1")
            change_background("none")
            os.rename("circos.svg", str(i) + ".svg")
            os.rename("circos.png", str(i) + ".png")
            if cogs:
                os.rename("temp/_prediction_deepnog.csv", "temp/" + str(i) + "_prediction_deepnog.csv")
            os.remove(file)
            i += 1
        merge.mergeImages(images, outFile = "circos.svg", align = circles_alignment, scale = scale, background_color = background_color)
        if legend or title != "":
            legendPosition = "top-right" if circles_alignment == "bottom" else "bottom-right"
            addText.addText(title, position = titlePos, inFile = "circos.svg", italic = italic, legend = legend, cogs_legend = cogs, legendPosition = legendPosition, cogs = full_cogs,
                            pCDS_color = CDS_positive, nCDS_color = CDS_negative, tRNA_color = tRNA, rRNA_color = rRNA, GC_content_color = GC_content)
            os.remove("circos.svg")
            os.rename("titled_circos.svg", "circos.svg")
    else:
        sizes, cogs_p, cogs_n = create_raw.base(gbk_file, "temp/", True, True, cogs, cogs, False, True, deepnog_bound)
        
        cogs_p = set(map(lambda x : "None" if x == None else x[0], cogs_p))
        cogs_n = set(map(lambda x : "None" if x == None else x[0], cogs_n))
        
        gbk2fna.gbkToFna(gbk_file, "temp/gbk_converted.fna")
        maxmins = GC_analysis.makeGC("temp/gbk_converted.fna", "temp/GC", window)
        createConf.create_conf(maxmins, GC_content, GC_skew, CDS_positive, CDS_negative, tRNA, skew_line_color, cogs, cogs_p=cogs_p, cogs_n=cogs_n)

        print("Drawing...")
        if which("circos") == None:
            print("Circos is not installed. please install for using this program.")
            raise(Exception)
        os.system("circos circos.conf >/dev/null 2>&1")
        os.system("circos -debug_group _all >/dev/null 2>&1")
        change_background(background_color)
        if legend or title != "":
            addText.addText(title, position = titlePos, inFile = "circos.svg", italic = italic, legend = legend, cogs_legend = cogs, legendPosition = "bottom-right", cogs = cogs_p.union(cogs_n),
            pCDS_color = CDS_positive, nCDS_color = CDS_negative, tRNA_color = tRNA, rRNA_color = rRNA, GC_content_color = GC_content)
            os.remove("circos.svg")
            os.rename("titled_circos.svg", "circos.svg")
    try:
        print("Converting to png...")
        from cairosvg import svg2png
        file = open("circos.svg")
        svg2png(bytestring = file.read(), write_to = output + ".png")
        file.close()
    except:
        print("There's been an error transforming to PNG, so the image might be incomplete or not exist. Please use the SVG instead. Did you install CairoSVG?")
    os.rename("circos.svg", output + ".svg")

    if not keep_temp_files:
        print("deleting temporary files")
        os.remove("circos.conf")
        for file in os.listdir("temp/"):
            os.remove("temp/" + file)
        for file in os.listdir("conf/"):
            os.remove("conf/" + file)
        # os.rmdir("temp") # don't know why I have permission denied for this.
    

def get_args():
    parser = ap.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str, help="Genbank file path", required=True)
    parser.add_argument("-o", "--output_file", type=str, help="Output image file path. Default: circos", default = "circos")
    parser.add_argument("-cu", "--cogs_unclassified", action='store_false', help="Do not clasify each coding sequence and draw them by color.", required = False)
    parser.add_argument("-b", "--deepnog_lower_bound", type=float, help="Lower bound for DeepNOG prediction certainty to be considered. Values in range [0,1] Default: 0", default = 0)
    parser.add_argument("-l", "--legend_not_included", action='store_false', help="Do not include color explanation.", required = False)
    parser.add_argument("-c", "-s", "--separate_circles", "--complete_genome", action='store_true', help="To draw each contig as a complete circle by itself.", required = False)
    parser.add_argument("-a", "--circles_alignment", type=str, choices=["center", "top", "bottom", "A", "<", "U"], help="When using --separate_circles, this defines the vertical alignment of every contig. Options: center, top, bottom, A (First on top), < (first to the left), U (Two on top, the rest below)", default = "auto")
    parser.add_argument("--scale", type=str, choices=["variable", "fixed", "sqrt"], help="When using --separate_circles, wether to use a different scale for tiny contigs, so to ensure visibility. Options: variable, fixed, sqrt", default = "sqrt")
    parser.add_argument("-k", "--keep_temporal_files", action='store_true', help="Don't delete files used for circos image generation, including protein categories prediction by Deepnog.", required = False)
    parser.add_argument("-w", "--window", "--step", type=int, help="base pair window for CG plotting. Default: 5000", default = 5000)
    title_group = parser.add_argument_group("title")
    title_group.add_argument("-t", "--title", type=str, help="Title of the image (strain name, or something like that). By default, it doesn't include title", default = "")
    title_group.add_argument("--title_position", type=str, choices=["center", "top", "bottom"], default = "center")
    title_group.add_argument("--italic_words", type=int, help="How many of the title's words should be written in italics. Default: 2", default = 2)
    color_group = parser.add_argument_group("colors")
    color_group.add_argument("-cs", "--color_scheme", "--color", type=str, choices=["auto", "blue", "purple", "soil", "greyscale", "velvet", "pastel"], help="Color scheme to use. Individual colors may be overriden wih other arguments. COGs coloring can't be changed.", default = 'auto')
    color_group.add_argument("-bc", "--background", "--background_color", type=str, help="Color for background. Default: none", default = 'none')
    color_group.add_argument("-pc", "--CDS_positive_color", type=str, help="Color for positive CDSs, in R, G, B format. Default: '180, 205, 222'", default = 'auto')
    color_group.add_argument("-nc", "--CDS_negative_color", type=str, help="Color for negative CDSs, in R, G, B format. Default: '53, 176, 42'", default = 'auto')
    color_group.add_argument("-tc", "--tRNA_color", type=str, help="Color for tRNAs, in R, G, B format. Default: '150, 5, 50'", default = 'auto')
    color_group.add_argument("-rc", "--rRNA_color", type=str, help="Color for rRNAs, in R, G, B format. Default: '150, 150, 50'", default = 'auto')
    color_group.add_argument("-cc", "--GC_content_color", type=str, help="Color for GC content, in R, G, B format. Default: '23, 0, 115'", default = "auto")
    color_group.add_argument("-sc", "--GC_skew_color", type=str, help="Color scheme for GC skew. For details on this, please read CIRCOS documentation. Default: 'eval(sprintf(\"rdbu-7-div-%%d\",remap_int(var(value),0,0,7,5)))'", default = 'auto')
    color_group.add_argument("-sl", "--GC_skew_line_color", type=str, help="Color for GC skew line. Default: black", default = 'auto')
    

    args = parser.parse_args()

    if args.output_file[-4:] == ".svg" or args.output_file[-4:] == ".png":
        args.output_file = args.output_file[:-4]
    if args.deepnog_lower_bound > 1 or args.deepnog_lower_bound < 0:
        print("DeepNOG lower bound must be between 0 and 1")
        raise Exception()

    return (args.input_file, args.output_file,
    args.cogs_unclassified, args.deepnog_lower_bound, args.legend_not_included, args.separate_circles, args.circles_alignment, args.scale, args.keep_temporal_files, args.window,
    args.title, args.title_position, args.italic_words, 
    args.color_scheme, args.background, args.GC_content_color, args.GC_skew_color, args.tRNA_color, args.rRNA_color, args.CDS_positive_color, args.CDS_negative_color, args.GC_skew_line_color)

if __name__ == "__main__":
    visualizeGenome(*get_args())
