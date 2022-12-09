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
            pCDS_colour = "180, 205, 222", nCDS_colour = "150, 200, 150", tRNA_colour = "150, 5, 50", rRNA_colour = "150, 150, 50", GC_content_colour = "23, 0, 115", size = "", font_colour = "0, 0, 0",
            tracks_explain = False):
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
            if tracks_explain:
                destiny.write('''
                    <text x="1500" y="535" font-size="30.0px" font-family="CMUBright-Roman" style="text-anchor:middle;" fill="{font_colour}">RNA (+)</text>
                    <text x="1500" y="570" font-size="30.0px" font-family="CMUBright-Roman" style="text-anchor:middle;" fill="{font_colour}">RNA (-)</text>
                    <text x="1500" y="495" font-size="50.0px" font-family="CMUBright-Roman" style="text-anchor:middle;" fill="{font_colour}">CDS (+)</text>
                    <text x="1500" y="630" font-size="50.0px" font-family="CMUBright-Roman" style="text-anchor:middle;" fill="{font_colour}">CDS (-)</text>
                    <text x="1500" y="810" font-size="50.0px" font-family="CMUBright-Roman" style="text-anchor:middle;" fill="{font_colour}">GC</text>
                    <text x="1500" y="840" font-size="35.0px" font-family="CMUBright-Roman" style="text-anchor:middle;" fill="{font_colour}">content</text>
                    <text x="1500" y="1050" font-size="50.0px" font-family="CMUBright-Roman" style="text-anchor:middle;" fill="{font_colour}">GC</text>
                    <text x="1500" y="1080" font-size="35.0px" font-family="CMUBright-Roman" style="text-anchor:middle;" fill="{font_colour}">skew</text>
                '''.format(font_colour=font_colour))
                if cogs_captions:
                    destiny.write('''
                        <text x="1500" y="440" font-size="50.0px" font-family="CMUBright-Roman" style="text-anchor:middle;" fill="{font_colour}">COG (+)</text>
                        <text x="1500" y="685" font-size="50.0px" font-family="CMUBright-Roman" style="text-anchor:middle;" fill="{font_colour}">COG (-)</text>
                    ''')
            if text != "":
                destiny.write(textElement)
            if size != "":
                destiny.write(sizeElement)
            if captions:

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
                    {"symbol": "X", "description":"[X] Mobilome: prophages, transposons", "colour": "83, 61, 145"}]
                COGs3 = [
                    {"symbol": "C", "description":"[C] Energy production and conversion", "colour": "181, 210, 123"},
                    {"symbol": "E", "description":"[E] Amino acid transport and metabolism", "colour": "131, 173, 41"},
                    {"symbol": "F", "description":"[F] Nucleotide transport and metabolism", "colour": "188, 184, 104"},
                    {"symbol": "G", "description":"[G] Carbohydrate transport and metabolism", "colour": "152, 143, 28"},
                    {"symbol": "H", "description":"[H] Coenzyme transport and metabolism", "colour": "237, 194, 138"},
                    {"symbol": "I", "description":"[I] Lipid transport and metabolism", "colour": "183, 130, 42"},
                    {"symbol": "P", "description":"[P] Inorganic ion transport and metabolism", "colour": "221, 137, 80"},
                    {"symbol": "Q", "description":"[Q] Secondary metabolities biosynthesis, transport and metabolism", "colour": "198, 95, 23"}]
                COGs4 = [
                    {"symbol": "R", "description":"[R] General function prediction only", "colour": "105, 105, 105"},
                    {"symbol": "S", "description":"[S] Function unknown", "colour": "153, 153, 153"},
                    {"symbol": "None", "description":"Unclassified", "colour": "234, 234, 234"}]
                
                if captionsPosition == "bottom":
                    legendElement = '<g style="transform: scale(3) translate(-2000px, -1700px)">'
                    X2 = 650.86
                elif captionsPosition == "top":
                    legendElement = '<g style="transform: scale(3) translate(-2000px, -300px)">'
                    X2 = 650.86
                else:
                    legendElement = '<g>'
                    X2 = 580.86

                allCogs = [ {"cogs": COGs1, "x": 233.95, "y": 44.68, "name": "Cellular Processes and Signaling"},
                            {"cogs": COGs2, "x": 233.95, "y": 162.9, "name": "Information Storage and Processing"},
                            {"cogs": COGs3, "x": X2, "y": 44.68, "name": "Metabolism"},
                            {"cogs": COGs4, "x": X2, "y": 162.9, "name": "Poorly Characterized"}
                ]

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
#                    '''elif captionsPosition == "top" or captionsPosition == "bottom":
#                        legendPiece = '''<rect x="{4}" y="{0}" width="30" height="30" fill="rgb({3})" stroke="{6}"/>
#                                <text x="{5}" y="{1}" font-size="38.0px" font-family="CMUBright-Roman" style="text-anchor:start;" fill="{6}">{2}</text>\n'''
#                        otherLegendPiece = '''<path d="M{4} {0} m0 20 h6 l4 -10 10 20 4-10 6 0" fill="{3}" stroke="{6}"/>
#                                <text x="{5}" y="{1}" font-size="40.0px" font-family="CMUBright-Roman" style="text-anchor:start;" fill="{6}">{2}</text>\n'''
#                        index = 0
#                        for COG in COGs1:
#                            
#                            if COG["symbol"] in cogs:
#                                legendElement += legendPiece.format((-900 if captionsPosition == "top" else 3000) + index*50, (-870 if captionsPosition == "top" else 3030) + index*50, COG["description"], COG["colour"],
#                                                 40, 80, font_colour)
#                                index += 1
#                            elif COG["symbol"] == "space":
#                                index += 1
#                        index = 0
#                        for COG in COGs2:
#                            
#                            if COG["symbol"] in cogs:
#                                legendElement += legendPiece.format((-900 if captionsPosition == "top" else 3000) + index*50, (-870 if captionsPosition == "top" else 3030) + index*50, COG["description"], COG["colour"],
#                                                 1340, 1380, font_colour)
#                                index += 1
#                            elif COG["symbol"] == "space":
#                                index += 1
#                        index = 0
#                        for feature in features:
#                            if feature["colour"] == "line":
#                                legendPiece = otherLegendPiece
#                                feature["colour"] = "none"
#                            if "eval" in feature["colour"]:
#                                feature["colour"] = "100,100,100"
#                            legendElement += legendPiece.format((-900 if captionsPosition == "top" else 3000) + index*50, (-870 if captionsPosition == "top" else 3030) + index*50, feature["description"], feature["colour"],
#                                             2640, 2680, font_colour)
#                            index += 1
#                    '''
                    else:
                        xb = 50 if "left" in captionsPosition else 2050
                        yb = 50 if "top" in captionsPosition else 2750
                        legendElement += '''<text font-size="20px" x="{x}" y="{y}">Genomic features</text>
                        <line stroke = "{font_colour}" x1="{xb}" y1="{yl}" x2="{xl}" y2="{yl}" />
                        <rect stroke = "{font_colour}" fill = "rgb({pCDS_color})" x="{xf}" y="{yf1}" width="18.66" height="18.66" />
                        <text x="{xft}" y="{yft1}" >positive CDS</text>
                        <rect stroke = "{font_colour}" fill = "rgb({nCDS_color})" x="{xf}" y="{yf2}" width="18.66" height="18.66" />
                        <text x="{xft}" y="{yft2}" >negative CDS</text>
                        <rect stroke = "{font_colour}" fill = "rgb({tRNA_color})" x="{xf}" y="{yf3}" width="18.66" height="18.66" />
                        <text x="{xft}" y="{yft3}">tRNA</text>
                        <rect stroke = "{font_colour}" fill = "rgb({rRNA_color})" x="{xf}" y="{yf4}" width="18.66" height="18.66" />
                        <text x="{xft}" y="{yft4}">rRNA</text>
                        <rect stroke = "{font_colour}" fill = "rgb({GC_color})" x="{xf}" y="{yf5}" width="18.66" height="18.66" />
                        <text x="{xft}" y="{yft5}">GC content</text>
                        <path stroke="{font_colour}" fill="none" d="M{xf6} {yf6} m0 10 h3 l2 -5 5 10 2-5 3 0" />
                        <text x="{xft}" y="{yft6}">GC skew</text>
                        '''.format( x = xb + 1.62,y = yb + 17.16,
                                    xb = xb, yl = yb + 22.52, xl = 167.56 + xb,
                                    xf = xb + 8.17, xft = xb + 35.5,
                                    yf1 = yb + 55.84, yf2 = yb + 82.73, yf3 = yb + 109.14, yf4 = yb + 136.03, yf5 = yb + 162.42, yf6 = yb + 192, xf6 = xb + 10,
                                    yft1 = yb + 70.54, yft2 = yb + 97.72, yft3 = yb + 124.89, yft4 = yb + 152.07, yft5 = yb + 179.25, yft6 = yb + 206.43,
                                    font_colour = font_colour,
                                    pCDS_color = pCDS_colour, nCDS_color = nCDS_colour, tRNA_color = tRNA_colour, rRNA_color = rRNA_colour, GC_color = GC_content_colour)

                        cogsAreDrawn = False
                        for functionalGroup in allCogs:
                            functionalGroupIsDrawn = False
                            i = -1
                            for cog in functionalGroup["cogs"]:
                                if cog["symbol"] in cogs:
                                    i += 1
                                    cogsAreDrawn = True
                                    functionalGroupIsDrawn = True
                                    x = functionalGroup["x"] + 30 + 87 * (i // 3)
                                    y = functionalGroup["y"] + 16 + 26 * (i % 3)
                                    legendElement += '''
                                        <rect stroke = "#13110C" fill="rgb({fill})" x="{x}" y="{y}" width="18.66" height="18.66"/>
                                        <text class="cls-2" x="{x2}" y="{y2}">[{name}] </text>'''.format(
                                        x = xb + x, y = yb + y, x2 = xb + x + 29, y2 = yb + y + 13, name = cog["symbol"], fill = cog["colour"]
                                    )
                            if functionalGroupIsDrawn:
                                legendElement += '''<text font-size="20px" x="{x}" y="{y}">{name}</text>'''.format(x = xb + functionalGroup["x"], y = yb + functionalGroup["y"], name = functionalGroup["name"])
                                legendElement += '''<line stroke = "#13110C" x1="{x1}" y1="{y}" x2="{x2}" y2="{y}" />'''.format(x1 = xb + functionalGroup["x"] - 0.36, x2 = xb + functionalGroup["x"] + 318, y = yb + functionalGroup["y"] + 5.4)
                        if cogsAreDrawn:
                            legendElement += '''<text font-size="20px" x="{x}" y="{y}">Cluster of Ortologues Groups (COGs)</text>'''.format(x = xb + 209.98, y = yb + 17.16)
                            legendElement += '''<line stroke = "#13110C" x1="{x1}" y1="{y}" x2="{x2}" y2="{y}" />'''.format(x1 = xb + 209.98 - 0.36, x2 = xb + 209.98 + 705, y = yb + 17.16 + 5.4)

                else: # if not cogs_legend
                    xb = 50 if "left" in captionsPosition else 2750
                    yb = 50 if "top" in captionsPosition else 2750
                    legendElement += '''<text font-size="20px" x="{x}" y="{y}">Genomic features</text>
                    <line stroke = "{font_colour}" x1="{xb}" y1="{yl}" x2="{xl}" y2="{yl}" />
                    <rect stroke = "{font_colour}" fill = "rgb({pCDS_color})" x="{xf}" y="{yf1}" width="18.66" height="18.66" />
                    <text x="{xft}" y="{yft1}" >positive CDS</text>
                    <rect stroke = "{font_colour}" fill = "rgb({nCDS_color})" x="{xf}" y="{yf2}" width="18.66" height="18.66" />
                    <text x="{xft}" y="{yft2}" >negative CDS</text>
                    <rect stroke = "{font_colour}" fill = "rgb({tRNA_color})" x="{xf}" y="{yf3}" width="18.66" height="18.66" />
                    <text x="{xft}" y="{yft3}">tRNA</text>
                    <rect stroke = "{font_colour}" fill = "rgb({rRNA_color})" x="{xf}" y="{yf4}" width="18.66" height="18.66" />
                    <text x="{xft}" y="{yft4}">rRNA</text>
                    <rect stroke = "{font_colour}" fill = "rgb({GC_color})" x="{xf}" y="{yf5}" width="18.66" height="18.66" />
                    <text x="{xft}" y="{yft5}">GC content</text>
                    <path stroke="{font_colour}" fill="none" d="M{xf6} {yf6} m0 10 h3 l2 -5 5 10 2-5 3 0" />
                    <text x="{xft}" y="{yft6}">GC skew</text>
                    '''.format( x = xb + 1.62,y = yb + 17.16,
                                xb = xb, yl = yb + 22.52, xl = 167.56 + xb,
                                xf = xb + 8.17, xft = xb + 35.5,
                                yf1 = yb + 55.84, yf2 = yb + 82.73, yf3 = yb + 109.14, yf4 = yb + 136.03, yf5 = yb + 162.42, yf6 = yb + 192, xf6 = xb + 10,
                                yft1 = yb + 70.54, yft2 = yb + 97.72, yft3 = yb + 124.89, yft4 = yb + 152.07, yft5 = yb + 179.25, yft6 = yb + 206.43,
                                font_colour = font_colour,
                                pCDS_color = pCDS_colour, nCDS_color = nCDS_colour, tRNA_color = tRNA_colour, rRNA_color = rRNA_colour, GC_color = GC_content_colour)

                legendElement += '</g>\n'
                destiny.write(legendElement)
        destiny.write(line)
    destiny.close()
    source.close()
