import argparse as ap
def create_conf(maxmins):
    file = open("circos.conf", "w")
    file.write('''karyotype = temp/output.kar
chromosomes_units = 100
chromosomes_display_default = yes

<<include conf/colors_fonts_patterns.conf>>

<<include conf/housekeeping.conf>>

<<include conf/highlight.conf>>

<<include conf/ticks.conf>>

# IMAGE
<image>
<<include image.conf>>
</image>

# IDEOGRAM
<ideogram>
<spacing>
default = 0.001r
</spacing>
radius           = 0.8r
thickness        = 40p
fill             = yes
stroke_color     = black
stroke_thickness = 2p
show_label       = no
label_radius     = dims(image,radius) - 250p
label_size       = 30p
label_parallel   = no
show_bands = yes
fill_bands = yes

</ideogram>

# PLOTS
<plots>

#Plot GC content
<plot>
type       = line
extend_bin = no
color      = black
fill_under = yes
thickness  = 1
file = temp/GC_GC_content.wig
r0   = 0.77r
r1   = 0.82r
min  = {GC_min}
max  = {GC_max}
<rules>
<rule>
condition = 1
fill_color = 224, 85, 113
</rule>
</rules>
</plot>

#Plot GC skew
<plot>
type       = line
extend_bin = no
color      = black
fill_under = no
thickness  = 1
file = temp/GC_GC_skew.wig
r0   = 0.68r
r1   = 0.75r
min  = {skew_min}
max  = {skew_max}
<rules>
<rule>
condition = 1
fill_color = 12, 6, 120
</rule>
</rules>
</plot>

# BACKGROUND
background                  = yes
background_stroke_color     = black
background_stroke_thickness = 1

</plots>
'''.format( GC_min = maxmins["min_GC_content"],
            GC_max = maxmins["max_GC_content"],
            skew_min = maxmins["min_skew"],
            skew_max = maxmins["max_skew"]))
    file.close()
    return