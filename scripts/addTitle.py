def addTitle(text, position = "center", inFile="VS4-2_circos.svg", outFile="default"):
    if(outFile == "default"):
        outFile = "titled_" + inFile
    textSize = 1500/len(text)
    if(position == "center"):
        verticalPosition = "1500"
    elif(position == "top"):
        verticalPosition = "150"
    elif(position == "bottom"):
        verticalPosition = "2900"
    textElement = '<text x="1500" y="{0}" font-size="{1}" font-family="CMUBright-Roman" text-anchor="middle">{2}</text>\n'.format(verticalPosition, textSize, text)
    # Can we calculate the best font size?
    source = open(inFile)
    destiny = open(outFile, "w")
    for line in source:
        if('</svg>' in line):
            destiny.write(textElement)
        destiny.write(line)
    destiny.close()
    source.close()
