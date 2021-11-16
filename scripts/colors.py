def parseColors(color_scheme = "auto", background_color = "none", GC_content = "auto", GC_skew ='auto', tRNA = 'auto', rRNA = 'auto', CDS_positive = 'auto', CDS_negative = 'auto', skew_line_color = '0, 0, 0'):
    if color_scheme == "blue" or color_scheme == "blues":
        GC_content = "14, 29, 130" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("rdbu-7-div-%d",remap_int(var(value),0,0,7,5)))' if GC_skew == "auto" else GC_skew
        tRNA = "99, 103, 186" if tRNA == "auto" else tRNA
        rRNA = "89, 123, 186" if rRNA == "auto" else rRNA
        CDS_positive = "191, 204, 217" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "171, 178, 217" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "163, 191, 217" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "purple" or color_scheme == "purples":
        GC_content = "57, 1, 120" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("puor-7-div-%d",remap_int(var(value),0,0,7,5)))' if GC_skew == "auto" else GC_skew
        tRNA = "108, 95, 156" if tRNA == "auto" else tRNA
        rRNA = "128, 85, 156" if rRNA == "auto" else rRNA
        CDS_positive = "196, 187, 237" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "151, 143, 191" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "163, 163, 217" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "soil":
        GC_content = "89, 60, 6" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("brbg-7-div-%d",remap_int(var(value),0,0,7,5)))' if GC_skew == "auto" else GC_skew
        tRNA = "0, 112, 11" if tRNA == "auto" else tRNA
        rRNA = "20, 112, 1" if rRNA == "auto" else rRNA
        CDS_positive = "185, 199, 186" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "175, 201, 178" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "89, 194, 115" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "greyscale" or color_scheme == "grayscale" or color_scheme == "grey" or color_scheme == "gray":
        GC_content = "87, 87, 87" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("rdgy-7-div-%d",remap_int(var(value),0,0,7,5)))' if GC_skew == "auto" else GC_skew
        tRNA = "115, 115, 115" if tRNA == "auto" else tRNA
        rRNA = "115, 115, 115" if rRNA == "auto" else rRNA
        CDS_positive = "209, 209, 209" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "184, 184, 184" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "171, 171, 171" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "velvet" or color_scheme == "pink" or color_scheme == "red":
        GC_content = "82, 1, 30" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("purd-7-seq-%d",remap_int(var(value),0,0,7,4)))' if GC_skew == "auto" else GC_skew
        tRNA = "130, 1, 49" if tRNA == "auto" else tRNA
        rRNA = "130, 21, 69" if rRNA == "auto" else rRNA
        CDS_positive = "217, 195, 205" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "196, 187, 193" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "230, 200, 211" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "pastel":
        GC_content = "196, 227, 255" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("pastel1-7-qual-%d",remap_int(var(value),0,0,4,3)))' if GC_skew == "auto" else GC_skew
        tRNA = "143, 143, 143" if tRNA == "auto" else tRNA
        rRNA = "143, 143, 143" if rRNA == "auto" else rRNA
        CDS_positive = "255, 222, 235" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "181, 255, 214" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "150, 150, 150" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "ocean" or color_scheme == "sea":
        GC_content = "24,138,141" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("%d,%d,127",remap_int(var(value),0,0,23,63),remap_int(var(value),0,0,87,159)))' if GC_skew == "auto" else GC_skew
        tRNA = "0, 83, 83" if tRNA == "auto" else tRNA
        rRNA = "0, 83, 83" if rRNA == "auto" else rRNA
        CDS_positive = "96,221,142" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "23,87,126" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "24,138,141" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "wood":
        GC_content = "150,90,49" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("%d,%d,%d",remap_int(var(value),0,0,120,163),remap_int(var(value),0,0,118,124),remap_int(var(value),0,0,37,81)))' if GC_skew == "auto" else GC_skew
        tRNA = "131,62,32" if tRNA == "auto" else tRNA
        rRNA = "62,131,32" if rRNA == "auto" else rRNA
        CDS_positive = "73,115,54" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "185,157,103" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "143,76,40" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "beach":
        GC_content = "231,207,133" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("%d,%d,%d",remap_int(var(value),0,0,93,210),remap_int(var(value),0,0,204,173),remap_int(var(value),0,0,167,128)))' if GC_skew == "auto" else GC_skew
        tRNA = "184,141,117" if tRNA == "auto" else tRNA
        rRNA = "117,167,184" if rRNA == "auto" else rRNA
        CDS_positive = "131,246,191" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "79,176,159" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "93,204,167" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "desert":
        GC_content = "225,168,111" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("%d,%d,%d",remap_int(var(value),0,0,237,193),remap_int(var(value),0,0,242,98),remap_int(var(value),0,0,154,2)))' if GC_skew == "auto" else GC_skew
        tRNA = "36,16,16" if tRNA == "auto" else tRNA
        rRNA = "180,165,150" if rRNA == "auto" else rRNA
        CDS_positive = "240,198,70" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "237,242,154" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "225,168,111" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "ice" or color_scheme == "frozen":
        GC_content = "186,242,239" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("%d,%d,%d",remap_int(var(value),0,0,37,162),remap_int(var(value),0,0,124,210),remap_int(var(value),0,0,163,223)))' if GC_skew == "auto" else GC_skew
        tRNA = "57,109,124" if tRNA == "auto" else tRNA
        rRNA = "97,129,144" if rRNA == "auto" else rRNA
        CDS_positive = "120,143,255" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "200,223,235" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "162,210,223" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "island":
        GC_content = "75,157,19" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("%d,%d,%d",remap_int(var(value),0,0,78,245),remap_int(var(value),0,0,247,229),remap_int(var(value),0,0,171,147)))' if GC_skew == "auto" else GC_skew
        tRNA = "75,157,19" if tRNA == "auto" else tRNA
        rRNA = "180,165,150" if rRNA == "auto" else rRNA
        CDS_positive = "254,155,227" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "154,215,173" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "155,187,89" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "forest":
        GC_content = "45,158,48" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("%d,%d,%d",remap_int(var(value),0,0,32,136),remap_int(var(value),0,0,103,208),remap_int(var(value),0,0,78,95)))' if GC_skew == "auto" else GC_skew
        tRNA = "109,137,75" if tRNA == "auto" else tRNA
        rRNA = "104,37,11" if rRNA == "auto" else rRNA
        CDS_positive = "136,208,95" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "105,184,128" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "45,158,48" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "toxic":
        GC_content = "137,4,177" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("%d,%d,%d",remap_int(var(value),0,0,106,11),remap_int(var(value),0,0,8,97),remap_int(var(value),0,0,136,11)))' if GC_skew == "auto" else GC_skew
        tRNA = "0,0,0" if tRNA == "auto" else tRNA
        rRNA = "50,50,50" if rRNA == "auto" else rRNA
        CDS_positive = "1,223,1" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "237,4,177" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "137,4,177" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "fire":
        GC_content = "216,0,0" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("%d,%d,0",remap_int(var(value),0,0,216,242),remap_int(var(value),0,0,0,85)))' if GC_skew == "auto" else GC_skew
        tRNA = "140,40,0" if tRNA == "auto" else tRNA
        rRNA = "161,0,0" if rRNA == "auto" else rRNA
        CDS_positive = "255,129,0" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "242,85,0" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "234,35,0" if skew_line_color == "auto" else skew_line_color
    elif color_scheme == "spring":
        GC_content = "249,150,174" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("%d,%d,%d",remap_int(var(value),0,0,217,107),remap_int(var(value),0,0,182,206),remap_int(var(value),0,0,253,238)))' if GC_skew == "auto" else GC_skew
        tRNA = "137,4,177" if tRNA == "auto" else tRNA
        rRNA = "237,164,207" if rRNA == "auto" else rRNA
        CDS_positive = "186,240,163" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "246,240,163" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "159,244,223" if skew_line_color == "auto" else skew_line_color
    else: # auto or neutral
        GC_content = "94, 120, 145" if GC_content == "auto" else GC_content
        GC_skew = 'eval(sprintf("bupu-7-seq-%d",remap_int(var(value),0,0,4,3)))' if GC_skew == "auto" else GC_skew
        tRNA = "67, 14, 110" if tRNA == "auto" else tRNA
        rRNA = "67, 110, 14" if rRNA == "auto" else rRNA
        CDS_positive = "186, 186, 186" if CDS_positive == "auto" else CDS_positive
        CDS_negative = "140, 140, 140" if CDS_negative == "auto" else CDS_negative
        skew_line_color = "171, 171, 171" if skew_line_color == "auto" else skew_line_color
    return color_scheme, background_color, GC_content, GC_skew, tRNA, rRNA, CDS_positive, CDS_negative, skew_line_color