import re

def addText(text, position = "center", inFile="circos.svg", outFile="default", italic=2,
            legend=True, legendPosition = "botom-right", cogs = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"}):
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
            textElement = '<text x="1500" y="{0}" font-size="{1}" font-family="CMUBright-Roman" text-anchor="middle">{2}</text>\n'.format(verticalPosition, textSize, text)
        else:
            textList = re.split(r"\s+", text)
            if len(textList) > italic:
                italicText = " ".join(textList[:italic])
                nonItalicText = " ".join(textList[italic:])
                center = (len(italicText) - len(nonItalicText)) * (5*float(textSize)/len(text)) + 1500
                print(center)
                textElement = '<text x="{4}" y="{0}" font-size="{1}" font-family="CMUBright-Roman" text-anchor="end" font-style="italic">{2}</text>\n<text x="{5}" y="{0}" font-size="{1}" font-family="CMUBright-Roman" text-anchor="start">  {3}</text>'.format(verticalPosition, textSize, italicText, nonItalicText, center-4, center+4)
            else:
                textElement = '<text x="1500" y="{0}" font-size="{1}" font-family="CMUBright-Roman" text-anchor="middle" font-style="italic">{2}</text>\n'.format(verticalPosition, textSize, text)

    source = open(inFile)
    destiny = open(outFile, "w")
    for line in source:
        if('</svg>' in line):
            if text != "":
                destiny.write(textElement)
            if legend:
                COGs1 = [{"symbol": "D", "description":"[D] Cell cycle control, adivision, chromosome partitioning", "color": "#697BB7"},
                {"symbol": "M", "description":"[M] Cell wall/membrane/envelope biogenesis", "color": "#304E9D"},
                {"symbol": "N", "description":"[N] Cell motility", "color": "#6A9AB5"},
                {"symbol": "O", "description":"[O] Post translational modification, protein turnover, chaperones", "color": "#1D82B0"},
                {"symbol": "T", "description":"[T] Signal transduction mechanism", "color": "#4CA9B5"},
                {"symbol": "U", "description":"[U] Intracellular trafficking, secretion and vesicular transport", "color": "#1C6875"},
                {"symbol": "V", "description":"[V] Defense mechanism", "color": "#6DBFA4"},
                {"symbol": "W", "description":"[W] Extracellular structures", "color": "#1A936F"},
                {"symbol": "Y", "description":"[Y] Nuclear structure", "color": "#4EB160"},
                {"symbol": "Z", "description":"[Z] Cytoskeleton", "color": "#1C7633"}]
                
                COGs2 = [{"symbol": "A", "description":"[A] RNA processing and modification", "color": "#E34949"},
                {"symbol": "B", "description":"[B] Chromatin structure and dynamics", "color": "#CD1B1B"},
                {"symbol": "J", "description":"[J] Translation, ribosomal structure and biogenesis", "color": "#AD5B9F"},
                {"symbol": "K", "description":"[K] Transcription", "color": "#A3378C"},
                {"symbol": "L", "description":"[L] Replication, recombination and repair", "color": "#8F76B4"},
                {"symbol": "X", "description":"[X] Mobilome: prophages, transposons", "color": "#533D91"},
                {"symbol": "space"},
                {"symbol": "C", "description":"[C] Energy production and conversion", "color": "#B5D27B"},
                {"symbol": "E", "description":"[E] Amino acid transport and metabolism", "color": "#83AD29"},
                {"symbol": "F", "description":"[F] Nucleotide transport and metabolism", "color": "#BCB868"},
                {"symbol": "G", "description":"[G] Carbohydrate transport and metabolism", "color": "#988F1C"},
                {"symbol": "H", "description":"[H] Coenzyme transport and metabolism", "color": "#EDC28A"},
                {"symbol": "I", "description":"[I] Lipid transport and metabolism", "color": "#B7822A"},
                {"symbol": "P", "description":"[P] Inorganic ion transport and metabolism", "color": "#DD8950"},
                {"symbol": "Q", "description":"[Q] Secondary metabolities biosynthesis, transport and metabolism", "color": "#C65F17"},
                {"symbol": "space"},
                {"symbol": "R", "description":"[R] General function prediction only", "color": "#696969"},
                {"symbol": "S", "description":"[S] Function unknown", "color": "#999999"},
                {"symbol": "None", "description":"Unclassified", "color": "#EAEAEA"}]

                legendElement = '<g>'
                index = 0
                legendPiece = '''<rect x="{4}" y="{0}" width="15" height="15" fill="{3}"/>
                                <text x="{5}" y="{1}" font-size="21.0px" font-family="CMUBright-Roman" style="text-anchor:start;">{2}</text>'''
                for COG in COGs1:
                    if COG["symbol"] in cogs:
                        if legendPosition == "top-right":
                            legendElement += legendPiece.format(50 + index*35, 68 + index*35, COG["description"], COG["color"], 1800, 1830)
                        elif legendPosition == "bottom-right":
                            legendElement += legendPiece.format(2730 + index*25, 2748 + index*25, COG["description"], COG["color"], 1800, 1830)
                        elif legendPosition == "top-left":
                            legendElement += legendPiece.format(50 + index*25, 68 + index*25, COG["description"], COG["color"], 650, 680)
                        else:
                            legendElement += legendPiece.format(2730 + index*25, 2748 + index*25, COG["description"], COG["color"], 650, 680)
                        index += 1
                index = 4
                for COG in COGs2:
                    if COG["symbol"] in cogs:
                        if legendPosition == "top-right":
                            legendElement += legendPiece.format(50 + index*35, 68 + index*35, COG["description"], COG["color"], 2400, 2430)
                        elif legendPosition == "bottom-right":
                            legendElement += legendPiece.format(2430 + index*25, 2448 + index*25, COG["description"], COG["color"], 2400, 2430)
                        elif legendPosition == "top-left":
                            legendElement += legendPiece.format(50 + index*25, 68 + index*25, COG["description"], COG["color"], 50, 80)
                        else:
                            legendElement += legendPiece.format(2430 + index*25, 2448 + index*25, COG["description"], COG["color"], 50, 80)
                        index += 1
                    elif COG["symbol"] == "space":
                        index += 1
                legendElement += '</g>'
                destiny.write(legendElement)
        destiny.write(line)
    destiny.close()
    source.close()
