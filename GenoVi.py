import argparse as ap
import scripts.create_raw as create_raw
import scripts.GC_analysis as GC_analysis
import os
import scripts.genbank2fna as gbk2fna
import scripts.createConf as createConf
import scripts.addText as addText
import scripts.mergeImages as merge
import glob
#from cairosvg import svg2png # The import is actually lower in the script and it is not always strictly necessary, so this should stay commented.

def visualizeGenome(gbk_file, output = "circos", 
                    cogs = True, legend = True, separate = False, circles_alignment = "center",
                    title = "", titlePos = "center", italic = 2,
                    GC_content = "23, 0, 115", GC_skew ='eval(sprintf("rdbu-7-div-%d",remap_int(var(value),0,0,7,5)))', tRNA = '150, 5, 50', CDS_positive = '180, 205, 222', CDS_negative = '53, 176, 42'):
    if not os.path.exists("temp"):
        os.mkdir("temp")
    legend = legend and not cogs
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
        for contig in contigs:
            file = "temp/" + str(i) + ".gbk"
            sizes, cogs_p, cogs_n = create_raw.base(file, "temp/", True, True, cogs, cogs, False)
            images.append({"size": sizes[0], "fileName": str(i) + ".svg"})
            gbk2fna.gbkToFna(file, "temp/gbk_converted.fna")
            maxmins = GC_analysis.makeGC("temp/gbk_converted.fna", "temp/GC")
            createConf.create_conf(maxmins, GC_content, GC_skew, CDS_positive, CDS_negative, tRNA, cogs, cogs_p, cogs_n)
    
            os.system("circos circos.conf")
            os.system("circos -debug_group _all")
            os.rename("circos.svg", str(i) + ".svg")
            os.rename("circos.png", str(i) + ".png")
            os.remove(file)
            i += 1
        merge.mergeImages(images, outFile = "circos.svg", align = circles_alignment)
    else:
        sizes, cogs_p, cogs_n = create_raw.base(gbk_file, "temp/", True, True, cogs, cogs, False)
        
        cogs_p = set(map(lambda x : x[0], cogs_p))
        cogs_n = set(map(lambda x : x[0], cogs_n))
        
        gbk2fna.gbkToFna(gbk_file, "temp/gbk_converted.fna")
        maxmins = GC_analysis.makeGC("temp/gbk_converted.fna", "temp/GC")
        createConf.create_conf(maxmins, GC_content, GC_skew, CDS_positive, CDS_negative, tRNA, cogs, cogs_p=cogs_p, cogs_n=cogs_n)

        os.system("circos circos.conf")
        os.system("circos -debug_group _all")
        if legend or title != "":
            addText.addText(title, position = titlePos, inFile = "circos.svg", italic = italic)
            os.remove("circos.svg")
            os.rename("titled_circos.svg", "circos.svg")
        os.rename("circos.png", output + ".png")
    try:
        from cairosvg import svg2png
        file = open("circos.svg")
        svg2png(bytestring = file.read(), write_to = output + ".png")
        file.close()
    except:
        print("We still haven't fully implemented png transformation, so the PNG image might be incomplete. Please use the SVG instead.")
    os.rename("circos.svg", output + ".svg")

    print("deleting temporary files")
    os.remove("temp/_bands.kar")
    os.remove("temp/_CDS_pos.txt")
    os.remove("temp/_CDS_neg.txt")
    os.remove("temp/_tRNA_pos.txt")
    os.remove("temp/_tRNA_neg.txt")
    os.remove("temp/gbk_converted.fna")
    os.remove("temp/GC_GC_skew.wig")
    os.remove("temp/GC_GC_content.wig")
    os.remove("circos.conf")
    
    for file in os.listdir("temp/"):
        os.remove("temp/" + file)
    # os.rmdir("temp") # don't know why I have permission denied for this.
    

def get_args():
    parser = ap.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str, help="Genbank file path", required=True)
    parser.add_argument("-o", "--output_file", type=str, help="Output image file path. Default: circos", default = "circos")
    parser.add_argument("-c", "--cogs_unclassified", action='store_false', help="Do not clasify each coding sequence and draw them by color.", required = False)
    parser.add_argument("-l", "--legend_not_included", action='store_false', help="Do not include color explanation.", required = False)
    parser.add_argument("-s", "--separate_circles", action='store_true', help="To draw each contig as a complete circle by itself.", required = False)
    parser.add_argument("-a", "--circles_alignment", type=str, choices=["center", "top", "bottom"], help="When using --separate_circles, this defines the vertical alignment of every contig. Options: center, top bottom", default = "center")
    parser.add_argument("-k", "--keep_temporal_files", action='store_true', help="Don't delete files used for circos image generation, including protein categories prediction by Deepnog.", required = False)
    title_group = parser.add_argument_group("title")
    title_group.add_argument("-t", "--title", type=str, help="Title of the image (strain name, or something like that). By default, it doesn't include title", default = "")
    title_group.add_argument("--title_position", type=str, choices=["center", "top", "bottom"], default = "center")
    title_group.add_argument("--italic_words", type=int, help="How many of the title's words should be written in italics. Default: 2", default = 2)
    color_group = parser.add_argument_group("colors")
    color_group.add_argument("-pc", "--CDS_positive_color", type=str, help="Color for positive CDSs, in R, G, B format. Default: '180, 205, 222'", default = '180, 205, 222')
    color_group.add_argument("-nc", "--CDS_negative_color", type=str, help="Color for negative CDSs, in R, G, B format. Default: '53, 176, 42'", default = '53, 176, 42')
    color_group.add_argument("-tc", "--tRNA_color", type=str, help="Color for tRNAs, in R, G, B format. Default: '150, 5, 50'", default = '150, 5, 50')
    color_group.add_argument("-cc", "--GC_content_color", type=str, help="Color for GC content, in R, G, B format. Default: '23, 0, 115'", default = "23, 0, 115")
    color_group.add_argument("-sc", "--GC_skew_color", type=str, help="Color scheme for GC skew. For details on this, please read CIRCOS documentation. Default: 'eval(sprintf(\"rdbu-7-div-%%d\",remap_int(var(value),0,0,7,5)))'", default = 'eval(sprintf("rdbu-7-div-%d",remap_int(var(value),0,0,7,5)))')
    

    args = parser.parse_args()

    if args.output_file[-4:] == ".svg":
        args.output_file[:-4]

    return (args.input_file, args.output_file,
    args.cogs_unclassified, args.legend_not_included, args.separate_circles, args.circles_alignment,
    args.title, args.title_position, args.italic_words, 
    args.GC_content_color, args.GC_skew_color, args.tRNA_color, args.CDS_positive_color, args.CDS_negative_color)

if __name__ == "__main__":
    visualizeGenome(*get_args())
