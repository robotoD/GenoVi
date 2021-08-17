import argparse as ap
import scripts.create_raw as create_raw
import scripts.GC_analysis as GC_analysis
import os
import scripts.genbank2fna as gbk2fna
import scripts.createConf as createConf
import scripts.addTitle as addTitle

def get_args():
    parser = ap.ArgumentParser()
    parser.add_argument("-i", "--input_file", type=str, help="Genbank file path", required=True)
    parser.add_argument("-o", "--output_file", type=str, help="Output image file path. Default: image", default = "image")
    parser.add_argument("-t", "--title", type=str, help="Title of the image (strain name, or something like that). Default: no title", default = "")
    parser.add_argument("--title_position", type=str, choices=["center", "top", "bottom"], default = "center")

    args = parser.parse_args()
        
    return args.input_file, args.output_file, args.title, args.title_position

if __name__ == "__main__":
    gbk_file, output, title, titlePos = get_args()[:]

    sizes = create_raw.create_kar(gbk_file, "temp/output.kar")
    create_raw.create_feature(gbk_file, "temp/cds_pos", "temp/cds_neg", sizes, "CDS")
    create_raw.create_feature(gbk_file, "temp/trna_pos", "temp/trna_neg", sizes, "tRNA")
    gbk2fna.gbkToFna(gbk_file, "temp/gbk_converted.fna")
    maxmins = GC_analysis.makeGC("temp/gbk_converted.fna", "temp/GC")
    createConf.create_conf(maxmins)
    os.system("circos circos.conf")
    os.system("circos -debug_group _all")
    if title != "":
        addTitle.addTitle(title, position = titlePos, inFile = "circos.svg")
        os.remove("circos.svg")
    
    print("deleting temporary files")
    os.remove("temp/output.kar")
    os.remove("temp/cds_pos")
    os.remove("temp/cds_neg")
    os.remove("temp/trna_pos")
    os.remove("temp/trna_neg")
    os.remove("temp/gbk_converted.fna")
    os.remove("temp/GC_GC_skew.wig")
    os.remove("temp/GC_GC_content.wig")
    os.remove("circos.conf")
    
    
    