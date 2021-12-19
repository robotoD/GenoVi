# GenoVi is a pipeline that generates circular maps for bacterial (complete or non-complete)
# genomes using Circos software. It also allows the user to annotate COG classifications
# through DeepNOG predictions.
# 
# GenoVi is under a BY-NC-SA Creative Commons License, Please cite. Cumsille et al., 2021
# You may remix, tweak, and build upon this work even for commercial purposes, as long as
# you credit this work and license your new creations under the identical terms.
# 
# Developed by Andres Cumsille, Andrea Rodriguez, Roberto E. Duran & Vicente Saona Urmeneta
# For any code related query, contact: andrea.rodriguezdelherbe@rdm.ox.ac.uk, vicente.saona@sansano.usm.cl.
#
# This script generates an image including all the others (Circos-generated SVGs)
# Input: List of dictionaries, with filename and image desired size (relative)
# Example: [{"fileName": "img1.svg", "size": 30000},
#           {"fileName": "img2.svg", "size": 10000}]


from math import sqrt


def mergeImages(images, outFile = "merged.svg", align = "auto", scale = "variable", background_color = "none", sort = False, captions_position = "normal"):
    totalWidth = 0
    extraElements = ""
    if sort:
        images.sort(key=lambda x:x["size"], reverse = True)
    for image in images:
        totalWidth += image["size"]
    if scale == "variable":
        totalVariableWidth = 0
        rectangleSize = 600 # Must be in range [0-3000]
        rectIsDrawn = False
        for i in range(len(images)):
            images[i]["scale"] = "fixed"
            if images[i]["size"] < totalWidth/10:
                
                for image in images[i:]:
                    image["scale"] = "variable"
                    totalVariableWidth += image["size"]
                    totalWidth -= image["size"]
                totalWidth += totalWidth/4
                break
    elif scale == "sqrt":
        totalWidth = 0
        for i in range(len(images)):
            images[i]["size"] =  sqrt((images[i]["size"] * images[0]["size"]))
            totalWidth += images[i]["size"]
                
    
    file = open(outFile, "w")
    file.write('''<?xml version="1.0" encoding="utf-8" standalone="no"?>
                <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">''')
    if captions_position == "right":
        file.write('''<svg width="4000px" height="3000px"''')
    elif captions_position == "left":
        file.write('''<svg width=viewBox="-1000 0 4000 3000"''')
    else:
        file.write('''<svg width="3000px" height="3000px"''')
    file.write('''version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">\n\n''')
    file.write('<rect x="-1000" y="0" width="5000px" height="3000px" style="fill:{};"/>'.format(background_color))
    currentX = 0
    beginGroup = '<g transform="translate({},{}) scale({})\n">'
    if align == "auto":
        if len(images) > 1 and images[1]["size"] > images[0]["size"]*0.75:
            align = "U"
        else:
            align = "A"
    imageIndex = 0
    for image in images:
        inFile = open(image["fileName"])
        for line in inFile:
            if("<svg " in line):
                break
        if align == "center":
            y= 3000 * (totalWidth - image["size"]) / (2 * totalWidth)
            file.write(beginGroup.format(currentX, y, float(image["size"])/totalWidth))
        elif align == "bottom":
            y= 3000 * (totalWidth - image["size"]) / (totalWidth)
            file.write(beginGroup.format(currentX, y, float(image["size"])/totalWidth))
        elif align == "top":
            y= 0
            file.write(beginGroup.format(currentX, y, float(image["size"])/totalWidth))
        elif align == "A":
            if imageIndex == 0:
                firstSize = float(image["size"])/totalWidth
                file.write(beginGroup.format(1500 - (3000*firstSize/2.0), 0, firstSize))
                currentX += 1500 - 1500*(1-firstSize) - 3000*float(image["size"])/totalWidth
            else:
                if scale == "variable" and image["scale"] == "variable":
                    if not rectIsDrawn:
                        rectangleSize = (2950 - firstSize*3000) * totalVariableWidth / max([x["size"] for x in images[imageIndex:]])
                        extraElements = '''<rect x="{0}" y="{1}" width="{2}" height="{3}" fill="none" stroke = "black"/>'''.format(currentX, 3000*firstSize, rectangleSize, 35 + max([x["size"] for x in images[imageIndex:]])*(rectangleSize*11/12)/totalVariableWidth)
                        extraElements += '<text x="{0}" y="{1}" font-size="30" font-family="CMUBright-Roman" text-anchor="left">scale: x{2}</text>\n'.format(currentX, 3000*firstSize + max([x["size"] for x in images[imageIndex:]])*(rectangleSize*11/12)/totalVariableWidth + 30, round(totalWidth/((3000/rectangleSize)*totalVariableWidth), 1))
                        rectIsDrawn = True
                    file.write(beginGroup.format(currentX, 3000*firstSize, float(image["size"])/((3000/rectangleSize)*totalVariableWidth)))
                    image["size"] = image["size"]*totalWidth/((3000/rectangleSize)*totalVariableWidth)
                else:
                    file.write(beginGroup.format(currentX, 3000*firstSize, float(image["size"])/totalWidth))
        elif align == "<":
            if currentX == 0:
                firstSize = float(image["size"])/totalWidth
                file.write(beginGroup.format(0, 1500 - (3000*firstSize/2.0), firstSize))
                currentX += 1500 - 1500*(1-firstSize) - 3000*float(image["size"])/totalWidth
            else:
                if scale == "variable" and image["scale"] == "variable":
                    if not rectIsDrawn:
                        extraElements = '<rect x="{0}" y="{1}" width="{2}" height="{3}" fill="none" stroke = "black"/>'.format(3000*firstSize, currentX, 35 + max([x["size"] for x in images[imageIndex:]])*(rectangleSize*11/12)/totalVariableWidth, rectangleSize)
                        extraElements += '<text x="{0}" y="{1}" font-size="30" font-family="CMUBright-Roman" text-anchor="left">scale: x{2}</text>\n'.format(3000*firstSize,currentX + (rectangleSize + 35), round(totalWidth/((3000/rectangleSize)*totalVariableWidth)))
                        rectIsDrawn = True
                    file.write(beginGroup.format(3000*firstSize, currentX, float(image["size"])/((3000/rectangleSize)*totalVariableWidth)))
                    image["size"] = image["size"]*totalWidth/((3000/rectangleSize)*totalVariableWidth)
                else:
                    file.write(beginGroup.format(3000*firstSize, currentX, float(image["size"])/totalWidth))
        elif align == "U":
            if imageIndex == 0:
                firstSize = float(image["size"])/totalWidth
                file.write(beginGroup.format(1500 - (3000*firstSize), 0, firstSize))
                currentX = -1 - 3000*float(image["size"])/totalWidth
            elif imageIndex == 1:
                secondSize = float(image["size"])/totalWidth
                file.write(beginGroup.format(1500, 0, secondSize))
                currentX = 1500 + 1500*(firstSize + secondSize - 1) - 3000*float(image["size"])/totalWidth
            else:
                if scale == "variable" and image["scale"] == "variable":
                    if not rectIsDrawn:
                        extraElements = '<rect x="{0}" y="{1}" width="{2}" height="{3}" fill="none" stroke = "black"/>'.format(currentX, 3000*max([firstSize, secondSize]), rectangleSize, 35 + max([x["size"] for x in images[imageIndex:]])*(rectangleSize*11/12)/totalVariableWidth)
                        extraElements += '<text x="{0}" y="{1}" font-size="30" font-family="CMUBright-Roman" text-anchor="left">scale: x{2}</text>\n'.format(currentX, 3000*max([firstSize, secondSize]) + max([x["size"] for x in images[imageIndex:]])*(rectangleSize*11/12)/totalVariableWidth + 30, round(totalWidth/((3000/rectangleSize)*totalVariableWidth)))
                        rectIsDrawn = True
                    file.write(beginGroup.format(currentX, 3000*max([firstSize, secondSize]), float(image["size"])/((3000/rectangleSize)*totalVariableWidth)))
                    image["size"] = image["size"]*totalWidth/((3000/rectangleSize)*totalVariableWidth)
                else:
                    file.write(beginGroup.format(currentX, 3000*max([firstSize, secondSize]), float(image["size"])/totalWidth))
        for line in inFile:
            if("</svg>" in line):
                break
            file.write(line)
        file.write('</g>\n')
        inFile.close()
        currentX += 3000*float(image["size"])/totalWidth
        imageIndex += 1
    file.write(extraElements)
    file.write('</svg>')
    file.close()
    
