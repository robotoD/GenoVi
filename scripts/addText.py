# GenoVi is a pipeline that generates circular maps for bacterial (complete or non-complete)
# genomes using Circos software. It also allows the user to annotate COG classifications
# through DeepNOG predictions.
# 
# GenoVi is under a BY-NC-SA Creative Commons License, Please cite. Cumsille et al., 2021
# You may remix, tweak, and build upon this work even for commercial purposes, as long as
# you credit this work and license your new creations under the identical terms.
# 
# Developed by Andres Cumsille, Andrea Rodriguez, Roberto E. Duran & Vicente Saona Urmeneta
# For any code related query, contact: andrea.rodriguezdelherbe@rdm.ox.ac.uk, vicente.saona@sansano.usm.cl


import re

__all__ = ['addText',
           ]

# Adds text elements to Circos-generated SVG image.
# May include title, contig size and colour legend.
def addText(text, position = "center", inFile="circos.svg", outFile="default", italic=2,
            captions=True, cogs_captions=True, captionsPosition = "botom-right", cogs = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"},
            pCDS_colour = "180, 205, 222", nCDS_colour = "150, 200, 150", tRNA_colour = "150, 5, 50", rRNA_colour = "150, 150, 50", GC_content_colour = "23, 0, 115", size = "", font_colour = "0, 0, 0"):
    if re.match("^\s*[012]?\d?\d\s*,\s*[012]?\d?\d\s*,\s*[012]?\d?\d\s*$", font_colour):
        font_colour = "rgb(" + font_colour + ")"
    if(outFile == "default"):
        outFile = "titled_" + inFile
    if text != "":
        textSize = 1500/len(text)
        if(position == "center"):
            verticalPosition = "1500"
        elif(position == "top"):
            verticalPosition = "150"
        elif(position == "bottom"):
            verticalPosition = "2900"
        if italic==0:
            textElement = '<text x="1500" y="{0}" font-size="{1}" font-family="CMUBright-Roman" text-anchor="middle" fill="{3}">{2}</text>\n'.format(verticalPosition, textSize, text, font_colour)
        else:
            textList = re.split(r"\s+", text)
            if len(textList) > italic:
                italicText = " ".join(textList[:italic])
                nonItalicText = " ".join(textList[italic:])
                textElement = '<text x="1500" y="{0}" font-size="{1}" font-family="CMUBright-Roman" text-anchor="middle" fill="{4}"><tspan font-style="italic">{2} </tspan>{3}</text>\n'.format(verticalPosition, textSize, italicText, nonItalicText, font_colour)
            else:
                textElement = '<text x="1500" y="{0}" font-size="{1}" font-family="CMUBright-Roman" text-anchor="middle" font-style="italic" fill = "{3}">{2}</text>\n'.format(verticalPosition, textSize, text, font_colour)
    if size != "":
        if size > 1000000:
            size = str(round(size / 1000000, 2)) + " Mb"
        elif size > 1000:
            size = str(round(size / 1000, 2)) + " kb"
        sizeElement = '<text x="1500" y="1560" font-size="64" font-family="CMUBright-Roman" text-anchor="middle" fill="{1}">{0}</text>\n'.format(size, font_colour)
    source = open(inFile)
    destiny = open(outFile, "w")
    background_has_been_amplified = False
    for line in source:
        if captionsPosition == "right" and '<svg width="3000px" height="3000px"' in line:
            destiny.write('<svg width="4000px" height="3000px" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">\n')
            continue
        elif captionsPosition == "left" and '<svg width="3000px" height="3000px"' in line:
            destiny.write('<svg viewBox="-1000 0 4000 3000" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">\n')
            continue
        elif captionsPosition == "bottom" and '<svg width="3000px" height="3000px"' in line:
            destiny.write('<svg viewBox="0 0 3000 4000" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">\n')
            continue
        elif captionsPosition == "top" and '<svg width="3000px" height="3000px"' in line:
            destiny.write('<svg viewBox="0 -1000 3000 4000" version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">\n')
            continue
        elif not background_has_been_amplified and (captionsPosition == "right" or captionsPosition == "left" or captionsPosition == "top" or captionsPosition == "up" or captionsPosition == "bottom" or captionsPosition == "down") and re.match('<rect x="-?1?0?0?0" y="0" width="[34]000px" height="3000px" style="fill:(.*);"/>',line):
            fill = re.match('<rect x="-?1?0?0?0" y="-?1?0?0?0" width="[34]000px" height="[34]000px" style="fill:(.*);"/>',line).groups()[0]
            destiny.write('<rect x="-1000" y="-1000" width="5000px" height="5000px" style="fill:{};"/>\n'.format(fill))
            background_has_been_amplified = True
            continue
        if('</svg>' in line):
            if text != "":
                destiny.write(textElement)
            if size != "":
                destiny.write(sizeElement)
            if captions:
                legendElement = '<g>'
                index = 0
                legendPiece = '''<rect x="{4}" y="{0}" width="15" height="15" fill="rgb({3})" stroke="{6}"/>
                                <text x="{5}" y="{1}" font-size="16.0px" font-family="CMUBright-Roman" style="text-anchor:start;" fill="{6}">{2}</text>\n'''
                otherLegendPiece = '''<path d="M{4} {0} m0 10 h3 l2 -5 5 10 2-5 3 0" fill="{3}" stroke="{6}"/>
                                <text x="{5}" y="{1}" font-size="16.0px" font-family="CMUBright-Roman" style="text-anchor:start;" fill="{6}">{2}</text>\n'''

                features = [{"description":"positive CDS", "colour": pCDS_colour},
                    {"description":"negative CDS", "colour": nCDS_colour},
                    {"description":"tRNA", "colour": tRNA_colour},
                    {"description":"rRNA", "colour": rRNA_colour},
                    {"description":"GC content", "colour": GC_content_colour},
                    {"description":"GC skew", "colour": "line"}]
                if cogs_captions:
                    COGs1 = [{"symbol": "D", "description":"[D] Cell cycle control, adivision, chromosome partitioning", "colour": "105, 123, 183"},
                    {"symbol": "M", "description":"[M] Cell wall/membrane/envelope biogenesis", "colour": "48, 78, 157"},
                    {"symbol": "N", "description":"[N] Cell motility", "colour": "106, 154, 181"},
                    {"symbol": "O", "description":"[O] Post translational modification, protein turnover, chaperones", "colour": "29, 130, 176"},
                    {"symbol": "T", "description":"[T] Signal transduction mechanism", "colour": "76, 169, 181"},
                    {"symbol": "U", "description":"[U] Intracellular trafficking, secretion and vesicular transport", "colour": "28, 104, 117"},
                    {"symbol": "V", "description":"[V] Defense mechanism", "colour": "109, 191, 164"},
                    {"symbol": "W", "description":"[W] Extracellular structures", "colour": "26, 147, 111"},
                    {"symbol": "Y", "description":"[Y] Nuclear structure", "colour": "78, 177, 96"},
                    {"symbol": "Z", "description":"[Z] Cytoskeleton", "colour": "28, 118, 51"}]
                    
                    COGs2 = [{"symbol": "A", "description":"[A] RNA processing and modification", "colour": "227, 73, 73"},
                    {"symbol": "B", "description":"[B] Chromatin structure and dynamics", "colour": "205, 27, 27"},
                    {"symbol": "J", "description":"[J] Translation, ribosomal structure and biogenesis", "colour": "173, 91, 159"},
                    {"symbol": "K", "description":"[K] Transcription", "colour": "163, 55, 140"},
                    {"symbol": "L", "description":"[L] Replication, recombination and repair", "colour": "143, 118, 180"},
                    {"symbol": "X", "description":"[X] Mobilome: prophages, transposons", "colour": "83, 61, 145"},
                    {"symbol": "space"},
                    {"symbol": "C", "description":"[C] Energy production and conversion", "colour": "181, 210, 123"},
                    {"symbol": "E", "description":"[E] Amino acid transport and metabolism", "colour": "131, 173, 41"},
                    {"symbol": "F", "description":"[F] Nucleotide transport and metabolism", "colour": "188, 184, 104"},
                    {"symbol": "G", "description":"[G] Carbohydrate transport and metabolism", "colour": "152, 143, 28"},
                    {"symbol": "H", "description":"[H] Coenzyme transport and metabolism", "colour": "237, 194, 138"},
                    {"symbol": "I", "description":"[I] Lipid transport and metabolism", "colour": "183, 130, 42"},
                    {"symbol": "P", "description":"[P] Inorganic ion transport and metabolism", "colour": "221, 137, 80"},
                    {"symbol": "Q", "description":"[Q] Secondary metabolities biosynthesis, transport and metabolism", "colour": "198, 95, 23"},
                    {"symbol": "space"},
                    {"symbol": "R", "description":"[R] General function prediction only", "colour": "105, 105, 105"},
                    {"symbol": "S", "description":"[S] Function unknown", "colour": "153, 153, 153"},
                    {"symbol": "None", "description":"Unclassified", "colour": "234, 234, 234"}]
                    
                    if captionsPosition == "right" or captionsPosition == "left":
                        legendPiece = '''<rect x="{4}" y="{0}" width="15" height="15" fill="rgb({3})" stroke="{6}"/>
                                <text x="{5}" y="{1}" font-size="25.0px" font-family="CMUBright-Roman" style="text-anchor:start;" fill="{6}">{2}</text>\n'''
                        otherLegendPiece = '''<path d="M{4} {0} m0 10 h3 l2 -5 5 10 2-5 3 0" fill="{3}" stroke="{6}"/>
                                <text x="{5}" y="{1}" font-size="25.0px" font-family="CMUBright-Roman" style="text-anchor:start;" fill="{6}">{2}</text>\n'''
                        index = 0
                        for COG in [*COGs1, *COGs2]:
                            
                            if COG["symbol"] in cogs:
                                legendElement += legendPiece.format(50 + index*35, 64 + index*35, COG["description"], COG["colour"],
                                                 3000 if captionsPosition == "right" else -900, 3030 if captionsPosition == "right" else -870, font_colour)
                                index += 1
                            elif COG["symbol"] == "space":
                                index += 1
                        index += 1
                        for feature in features:
                            if feature["colour"] == "line":
                                legendPiece = otherLegendPiece
                                feature["colour"] = "none"
                            if "eval" in feature["colour"]:
                                feature["colour"] = "100,100,100"
                            legendElement += legendPiece.format(50 + index*35, 64 + index*35, feature["description"], feature["colour"],
                                             3000 if captionsPosition == "right" else -900, 3030 if captionsPosition == "right" else -870, font_colour)
                            index += 1
                    elif captionsPosition == "top" or captionsPosition == "bottom":
                        legendPiece = '''<rect x="{4}" y="{0}" width="30" height="30" fill="rgb({3})" stroke="{6}"/>
                                <text x="{5}" y="{1}" font-size="38.0px" font-family="CMUBright-Roman" style="text-anchor:start;" fill="{6}">{2}</text>\n'''
                        otherLegendPiece = '''<path d="M{4} {0} m0 20 h6 l4 -10 10 20 4-10 6 0" fill="{3}" stroke="{6}"/>
                                <text x="{5}" y="{1}" font-size="40.0px" font-family="CMUBright-Roman" style="text-anchor:start;" fill="{6}">{2}</text>\n'''
                        index = 0
                        for COG in COGs1:
                            
                            if COG["symbol"] in cogs:
                                legendElement += legendPiece.format((-900 if captionsPosition == "top" else 3000) + index*50, (-870 if captionsPosition == "top" else 3030) + index*50, COG["description"], COG["colour"],
                                                 40, 80, font_colour)
                                index += 1
                            elif COG["symbol"] == "space":
                                index += 1
                        index = 0
                        for COG in COGs2:
                            
                            if COG["symbol"] in cogs:
                                legendElement += legendPiece.format((-900 if captionsPosition == "top" else 3000) + index*50, (-870 if captionsPosition == "top" else 3030) + index*50, COG["description"], COG["colour"],
                                                 1340, 1380, font_colour)
                                index += 1
                            elif COG["symbol"] == "space":
                                index += 1
                        index = 0
                        for feature in features:
                            if feature["colour"] == "line":
                                legendPiece = otherLegendPiece
                                feature["colour"] = "none"
                            if "eval" in feature["colour"]:
                                feature["colour"] = "100,100,100"
                            legendElement += legendPiece.format((-900 if captionsPosition == "top" else 3000) + index*50, (-870 if captionsPosition == "top" else 3030) + index*50, feature["description"], feature["colour"],
                                             2640, 2680, font_colour)
                            index += 1
                    else:
                        for COG in COGs1[::-1]:
                            if COG["symbol"] in cogs:
                                if captionsPosition == "top-right":
                                    legendElement += legendPiece.format(50 + index*25, 62 + index*25, COG["description"], COG["colour"], 1800, 1830, font_colour)
                                elif captionsPosition == "bottom-right":
                                    legendElement += legendPiece.format(2980 - index*25, 2992 - index*25, COG["description"], COG["colour"], 1800, 1830, font_colour)
                                elif captionsPosition == "top-left":
                                    legendElement += legendPiece.format(50 + index*25, 62 + index*25, COG["description"], COG["colour"], 650, 680, font_colour)
                                else:
                                    legendElement += legendPiece.format(2980 - index*25, 2992 - index*25, COG["description"], COG["colour"], 650, 680, font_colour)
                                index += 1
                        index = 0
                        for COG in COGs2[::-1]:
                            if COG["symbol"] in cogs:
                                if captionsPosition == "top-right":
                                    legendElement += legendPiece.format(50 + index*25, 62 + index*25, COG["description"], COG["colour"], 2400, 2430, font_colour)
                                elif captionsPosition == "bottom-right":
                                    legendElement += legendPiece.format(2980 - index*25, 2992 - index*25, COG["description"], COG["colour"], 2400, 2430, font_colour)
                                elif captionsPosition == "top-left":
                                    legendElement += legendPiece.format(50 + index*25, 62 + index*25, COG["description"], COG["colour"], 50, 80, font_colour)
                                else:
                                    legendElement += legendPiece.format(2980 - index*25, 2992 - index*25, COG["description"], COG["colour"], 50, 80, font_colour)
                                index += 1
                            elif COG["symbol"] == "space":
                                index += 1
                        index = 0
                        for feature in features:
                            if feature["colour"] == "line":
                                legendPiece = otherLegendPiece
                                feature["colour"] = "none"
                            if "eval" in feature["colour"]:
                                feature["colour"] = "100,100,100"
                            if captionsPosition == "top-right":
                                legendElement += legendPiece.format(50 + index*25, 62 + index*25, feature["description"], feature["colour"], 1600, 1630, font_colour)
                            elif captionsPosition == "bottom-right":
                                legendElement += legendPiece.format(2855 + index*25, 2862 + index*25, feature["description"], feature["colour"], 1600, 1630, font_colour)
                            elif captionsPosition == "top-left":
                                legendElement += legendPiece.format(50 + index*25, 62 + index*25, feature["description"], feature["colour"], 450, 480, font_colour)
                            else:
                                legendElement += legendPiece.format(2855 + index*25, 2862 + index*25, feature["description"], feature["colour"], 450, 480, font_colour)
                            index += 1
                else: # if not cogs_legend
                    index = 0
                    if captionsPosition == "right":
                        legendPiece = '''<rect x="{4}" y="{0}" width="15" height="15" fill="rgb({3})" stroke="{6}"/>
                                <text x="{5}" y="{1}" font-size="25.0px" font-family="CMUBright-Roman" style="text-anchor:start;" fill="{6}">{2}</text>\n'''
                        otherLegendPiece = '''<path d="M{4} {0} m0 10 h3 l2 -5 5 10 2-5 3 0" fill="{3}" stroke="{6}"/>
                                <text x="{5}" y="{1}" font-size="25.0px" font-family="CMUBright-Roman" style="text-anchor:start;" fill="{6}">{2}</text>\n'''
                        for feature in features:
                            if feature["colour"] == "line":
                                legendPiece = otherLegendPiece
                                feature["colour"] = "none"
                            if "eval" in feature["colour"]:
                                feature["colour"] = "100,100,100"
                            legendElement += legendPiece.format(50 + index*35, 64 + index*35, feature["description"], feature["colour"], 3000, 3030, font_colour)
                            index += 1
                    elif captionsPosition == "left":
                        legendPiece = '''<rect x="{4}" y="{0}" width="15" height="15" fill="rgb({3})" stroke="{6}"/>
                                <text x="{5}" y="{1}" font-size="25.0px" font-family="CMUBright-Roman" style="text-anchor:start;" fill="{6}">{2}</text>\n'''
                        otherLegendPiece = '''<path d="M{4} {0} m0 10 h3 l2 -5 5 10 2-5 3 0" fill="{3}" stroke="{6}"/>
                                <text x="{5}" y="{1}" font-size="25.0px" font-family="CMUBright-Roman" style="text-anchor:start;" fill="{6}">{2}</text>\n'''
                        for feature in features:
                            if feature["colour"] == "line":
                                legendPiece = otherLegendPiece
                                feature["colour"] = "none"
                            if "eval" in feature["colour"]:
                                feature["colour"] = "100,100,100"
                            legendElement += legendPiece.format(50 + index*35, 64 + index*35, feature["description"], feature["colour"], -900, -870, font_colour)
                            index += 1
                    else:
                        for feature in features:
                            if feature["colour"] == "line":
                                legendPiece = otherLegendPiece
                                feature["colour"] = "none"
                            if "eval" in feature["colour"]:
                                feature["colour"] = "100,100,100"
                            if captionsPosition == "top-right":
                                legendElement += legendPiece.format(50 + index*25, 68 + index*25, feature["description"], feature["colour"], 2700, 2730, font_colour)
                            elif captionsPosition == "bottom-right":
                                legendElement += legendPiece.format(2830 + index*25, 2848 + index*25, feature["description"], feature["colour"], 2700, 2730, font_colour)
                            elif captionsPosition == "top-left":
                                legendElement += legendPiece.format(50 + index*25, 68 + index*25, feature["description"], feature["colour"], 50, 80, font_colour)
                            else:
                                legendElement += legendPiece.format(2830 + index*25, 2848 + index*25, feature["description"], feature["colour"], 50, 80, font_colour)
                            index += 1
                legendElement += '</g>\n'
                destiny.write(legendElement)
        destiny.write(line)
    destiny.close()
    source.close()
