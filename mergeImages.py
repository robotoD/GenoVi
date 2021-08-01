# Generates an image including all the others
# Input: List of dictionaries, with file name and image desired size (relative)
# Example: [{"fileName": "img1.svg", "size": 30000},
#           {"fileName": "img2.svg", "size": 10000}]
def mergeImages(images, outFile = "merged.svg"):
    totalWidth = 0
    images.sort(key=lambda x:x["size"], reverse = True)
    for image in images:
        totalWidth += image["size"]
    
    file = open(outFile, "w")
    file.write('''<?xml version="1.0" encoding="utf-8" standalone="no"?>
                <!DOCTYPE svg PUBLIC "-//W3C//DTD SVG 1.1//EN" "http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd">
                <svg width="3000px" height="3000px"
                version="1.1" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink">\n\n''')
    currentX = 0
    for image in images:
        inFile = open(image["fileName"])
        for line in inFile:
            if("<svg " in line):
                break
        file.write('<g transform="translate({},0) scale({})\n">'.format(currentX, float(image["size"])/totalWidth))
        for line in inFile:
            if("</svg>" in line):
                break
            file.write(line)
        file.write('</g>\n')
        inFile.close()
        currentX += 3000*float(image["size"])/totalWidth
    file.write('</svg>')
    file.close()
    
