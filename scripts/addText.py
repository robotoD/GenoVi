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
                COGs = ["1- COG", "A- COG", "A- COG", "A- COG", "5- COG", "A- COG", "A- COG", "A- COG", "9- COG",
                "10- COG", "A- COG", "A- COG", "A- COG", "A- COG", "15- COG", "16- COG", "17- COG", "18- COG", "19- COG", "20- COG"]
                colors = ["#FF0000", "#00FF00", "#0000FF", "#FFFF00", "#FF00FF", "#00FFFF", "#FFAA00", "#FF00AA", "#FFAAAA", "#AAFF00",
                "#AA00FF", "#AAFFAA", "#00FFAA", "#AA00FF", "#00AAFF", "#AAAAFF", "#AAAAAA", "#AA0000", "#00AA00", "#0000AA"]
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
