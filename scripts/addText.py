import re

def addText(text, position = "center", inFile="circos.svg", outFile="default", italic=2, legend=True):
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
                COGs = ["[D] Cell cycle control, adivision, chromosome partitioning",
                "[M] Cell wall/membrane/envelope biogenesis",
                "[N] Cell motility",
                "[O] Post translational mmodification, protein turnover, chaperones",
                "[T] Signal transduction mechanism",
                "[U] Intracellular trafficking, secretion and vesicular transport",
                "[V] Defense mechanism",
                "[W] Extracellular structures",
                "[Y] Nuclear structure",
                "[Z] Cytoskeleton",
                
                "[A] RNA processing and modification",
                "[B] Chromatin structure and dynamics",
                "[J] Translation, ribosomal structure and biogenesis",
                "[K] Transcription",
                "[L] Replication, recombination and repair",
                "[X] Mobilome: prophages, transposons",
                
                "[C] Energy production and conversion",
                "[E] Amino acid transport and metabolism",
                "[F] Nucleotide transport and metabolism",
                "[G] Carbohydrate transport and metabolism",
                "[H] Coenzyme transport and metabolism",
                "[I] Lipid transport and metabolism",
                "[P] Inorganic ion transport and metabolism",
                "[Q] Secondary metabolities biosynthesis, transport and metabolism",
                
                "[R] General function prediction only",
                "[S] Function unknown",
                "Unclassified"]
                colors = ["#697BB7", "#304E9D", "#6A9AB5", "#1D82B0", "#4CA9B5", "#1C6875", "#6DBFA4", "#1A936F", "#4EB160", "#1C7633",
                "#E34949", "#CD1B1B", "#AD5B9F", "#A3378C", "#8F76B4", "#533D91",
                "#B5D27B", "#83AD29", "#BCB868", "#988F1C", "#EDC28A", "#B7822A", "#DD8950", "#C65F17",
                "#696969", "#999999", "EAEAEA"]
                legendElement = '<g>'
                index = 0
                for i in range(len(COGs)):
                    legendElement += '''<rect x="2500" y="{0}" width="20" height="20" fill="{3}"/>
                                        <text x="2530" y="{1}" font-size="26.0px" font-family="CMUBright-Roman" style="text-anchor:start;">{2}</text>'''.format(2400 + index*30, 2418 + index*30, COGs[i], colors[i])
                    index += 1
                legendElement += '''</g>'''
                destiny.write(legendElement)
        destiny.write(line)
    destiny.close()
    source.close()
