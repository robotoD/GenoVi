import re

def addTitle(text, position = "center", inFile="circos.svg", outFile="default", italic=2):
    if(outFile == "default"):
        outFile = "titled_" + inFile
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
    # Can we calculate the best font size?
    source = open(inFile)
    destiny = open(outFile, "w")
    for line in source:
        if('</svg>' in line):
            destiny.write(textElement)
        destiny.write(line)
    destiny.close()
    source.close()

