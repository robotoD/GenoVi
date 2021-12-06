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


import argparse as ap
import os


# Function that writes base CIRCOS configuration files
# circos.conf, highlight.conf, colors_fonts_patterns.conf, housekeeping.conf, image.conf, and ticks.conf.
def create_conf(maxmins,
                font_color = "0, 0, 0",
                GC_content_color = "23, 0, 115",
                GC_skew_color = 'eval(sprintf("rdbu-7-div-%d",remap_int(var(value),0,0,7,5)))',
                CDS_positive_color = '180, 205, 222',
                CDS_negative_color = '53, 176, 42',
                tRNA_color = '150, 5, 50',
                rRNA_color = '150, 150, 50',
                GC_skew_line_color = '0, 0, 0',
                background_color = "transparent",
                cogs = True,
                cogs_p = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "None"},
                cogs_n = {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z", "None"}):
    file = open("circos.conf", "w")
    file.write('''karyotype = temp/_bands.kar
chromosomes_units = 100
chromosomes_display_default = yes

<<include conf/colors_fonts_patterns.conf>>

<<include conf/housekeeping.conf>>

<<include conf/highlight.conf>>

<<include conf/ticks.conf>>

# IMAGE
<image>
<<include conf/image.conf>>
</image>

# IDEOGRAM
<ideogram>
<spacing>
default = 0.001r
</spacing>
radius           = 0.8r
thickness        = 40p
fill             = yes
stroke_color     = {font_color}
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
color      = {font_color}
fill_under = yes
thickness  = 1
file = temp/GC_GC_content.wig
r0   = 0.62r
r1   = 0.7r
min  = {GC_min}
max  = {GC_max}
<rules>
<rule>
condition = 1
fill_color = {content_fill}
</rule>
</rules>
</plot>

#Plot GC skew
<plot>
type       = line
extend_bin = no
color      = {skew_line}
fill_under = no
thickness  = 0.5
file = temp/GC_GC_skew.wig
r0   = 0.3r
r1   = 0.6r
min  = {skew_min}
max  = {skew_max}
<rules>
<rule>
condition = 1
fill_color = {skew_fill}
</rule>
</rules>
</plot>

# BACKGROUND
background                  = yes
background_stroke_color     = {font_color}
background_stroke_thickness = 1

</plots>
'''.format( GC_min = maxmins["min_GC_content"],
            GC_max = maxmins["max_GC_content"],
            skew_min = maxmins["min_skew"],
            skew_max = maxmins["max_skew"],
            content_fill = GC_content_color,
            skew_fill = GC_skew_color,
            skew_line = GC_skew_line_color,
            font_color = font_color))
    file.close()
    if not os.path.exists("conf"):
        os.mkdir("conf")
    
    file = open("conf/highlight.conf", "w")
    file.write('''fill_color       = grey
color            = {font_color}
stroke_thickness = 1
r1               = 0.5r
r0               = 1r

<highlights>

<highlight>

#Positive band BGCs
init_counter = highlight:1
file = temp/_CDS_pos.txt
fill_color = {CDS_positive}
r1 = 0.85r
r0 = 0.89r
</highlight>

#Negative
<highlight>

file = temp/_CDS_neg.txt
fill_color = {CDS_negative}
r1 = 0.8r
r0 = 0.84r
</highlight>

#Positive band tRNAs
<highlight>
file = temp/_tRNA_pos.txt
fill_color = {tRNA_positive}
r1 = 0.85r
r0 = 0.89r
</highlight>

#Negative
<highlight>
file = temp/_tRNA_neg.txt
fill_color = {tRNA_negative}
r1 = 0.8r
r0 = 0.84r
</highlight>

#Positive band rRNAs
<highlight>
file = temp/_rRNA_pos.txt
fill_color = {rRNA_positive}
r1 = 0.85r
r0 = 0.89r
</highlight>

#Negative
<highlight>
file = temp/_rRNA_neg.txt
fill_color = {rRNA_negative}
r1 = 0.8r
r0 = 0.84r
</highlight>

'''.format(CDS_positive = CDS_positive_color,
           CDS_negative = CDS_negative_color,
           tRNA_positive = tRNA_color,
           tRNA_negative = tRNA_color,
           rRNA_positive = rRNA_color,
           rRNA_negative = rRNA_color,
           font_color = font_color))
    if cogs:
        for COG in [{"name": "D", "color": "99, 123, 183"},
                    {"name": "M", "color": "38, 89, 168"},
                    {"name": "N", "color": "100, 151, 176"},
                    {"name": "O", "color": "32, 132, 174"},
                    {"name": "T", "color": "66, 169, 179"},
                    {"name": "U", "color": "18, 105, 116"},
                    {"name": "V", "color": "97, 195, 166"},
                    {"name": "W", "color": "37, 151, 117"},
                    {"name": "Y", "color": "70, 185, 109"},
                    {"name": "Z", "color": "21, 121, 60"},

                    {"name": "A", "color": "231, 82, 81"},
                    {"name": "B", "color": "208, 36, 41"},
                    {"name": "J", "color": "166, 96, 167"},
                    {"name": "K", "color": "157, 63, 151"},
                    {"name": "L", "color": "133, 116, 181"},
                    {"name": "X", "color": "79, 72, 158"},
                    
                    {"name": "C", "color": "170, 211, 130"},
                    {"name": "E", "color": "125, 176, 64"},
                    {"name": "F", "color": "178, 179, 109"},
                    {"name": "G", "color": "143, 138, 47"},
                    {"name": "H", "color": "235, 188, 134"},
                    {"name": "I", "color": "175, 126, 53"},
                    {"name": "P", "color": "219, 136, 86"},
                    {"name": "Q", "color": "196, 100, 38"},
                    
                    {"name": "R", "color": "98, 98, 98"},
                    {"name": "S", "color": "144, 144, 144"},
                    {"name": "None", "color": "234, 234, 234"},
                    ]:
            if COG["name"] in cogs_p:
                file.write('''
#positive {name}
<highlight>

file = temp/_CDS_pos_{name}.txt
fill_color = {color}
r1 = 0.9r
r0 = 0.94r
</highlight>
'''.format(name = COG["name"], color = COG["color"]))
            if COG["name"] in cogs_n:
                file.write('''
#negative {name}
<highlight>

file = temp/_CDS_neg_{name}.txt
fill_color = {color}
r1 = 0.75r
r0 = 0.79r
</highlight>
'''.format(name = COG["name"], color = COG["color"]))
    file.write("</highlights>")
    file.close()

    file = open("conf/colors_fonts_patterns.conf", "w")
    file.write('''<colors>
<<include etc/colors.conf>>
</colors>

<fonts>
<<include etc/fonts.conf>>
</fonts>

<patterns>
<<include etc/patterns.conf>>
</patterns>
''')
    file.close()

    file = open("conf/housekeeping.conf", "w")
    file.write('''anglestep       = 0.5
minslicestep    = 10
beziersamples   = 40 # bezier curves are drawn piece-wise
                     # linear, with this many samples
debug           = no
warnings        = no
imagemap        = no
paranoid        = yes

units_ok        = bupr
units_nounit    = n

# \t  tab
# \s  any whitespace
file_delim = \s
# collapse adjacent whitespace 
# e.g. two spaces are treated as one, not as a missing field
file_delim_collapse = yes

# Record delimiter for parameter values that are lists, such as
# hs1:0.25;hs2:0.10. By default, both ; and , are accepted
#
# e.g. hs1:0.25,hs2:0.10
#      hs1:0.25;hs2:0.10
list_record_delim = \s*[;,]\s*
# Field delimiter specifies the assignment operator, e.g. 
list_field_delim  = \s*[:=]\s*

# Record delimiter for options in data files.
options_record_delim = [,;]
options_field_delim  = =

# Rule fields and other parameters accept var(VARIABLE) syntax
# to reference parameters of data points. By default, if
# VARIABLE does not exist, Circos quits with an error, unless
# the skip parameter below is set.
# 
# This feature is useful when you have data that don't always
# have the same options. For example,
# 
# chr1 10 20 a=10
# chr1 50 60 b=10
skip_missing_expression_vars = no

# In old versions, data point parameters were referenced using _NAME_
# syntax. This has been replaced with var(NAME). The _NAME_ syntax is
# deprecated (for example, it will break when dealing with fields like
# gene_a_1). If you must use it, set the parameter below.

legacy_underline_expression_syntax = no

# Magnification factor for text in SVG files.
svg_font_scale = 1.3

# super/sub script baseline shift (%)
sup_baseline_shift = 40
sub_baseline_shift = -40

# super/sub script font size (%)
sup_fontsize = 90
sub_fontsize = 90

# default font - pick one of the keys from <fonts> block
default_font   = default
# default font name is used for SVG files for cases where
# the font definition does not include a name
# see etc/fonts.conf for details
default_font_name  = Arial
default_font_color = {font_color}

# default color for cases when color is not specified
default_color  = {font_color}

<guides>
thickness      = 1
size           = 5
type           = outline
<object>
all            = no
ideogram       = no
ideogram_label = no
</object>
<color>
default = lblue
text    = red
</color>
</guides>

# Receive debug messages about actions
# 
# For a full list of available debug groups, see 
#   
#  http://www.circos.ca/documentation/tutorials/configuration/debugging/
#
# Common ones are:
#
# summary   - top level indications of what's happening (on by default)
# chrfilter - ideogram filtering (parsing 'chromosomes' parameter)
# conf      - configuration file
# io        - file I/O
# ideogram  - ideogram reporting
# layers    - specific plot z-layers
# rules     - dynamic rules
# color     - color allocation and definition
# timers    - some benchmark timings
# cache     - caches
# _all      - turn on all groups

debug_group = summary,output

# run length duration after which timing report is automatically
# generated at the end of the run
debug_auto_timer_report = 30

debug_word_separator = " "
debug_undef_text     = _undef_
debug_empty_text     = _emptylist_

# parameters passed to functions can be
# validated to check consistency. turn this
# off to speed things up
debug_validate       = yes

# Reformat numbers in debug output for consistency.
# If you have a lot of debug output (e.g. -debug_group _all)
# this will slow things considerably.
debug_output_tidy    = no

# pixel sub-sampling for text tracks
text_pixel_subsampling = 1
# array or span - use 'span' if applying snuggle refinement
text_snuggle_method    = array

# restrict names of parameters?
# if 'yes' then only certain named parameters are allowed within
# blocks and option fields for data
#
# if 'no' then you can define parameters with any name what-so-ever,
# useful if you wish to define states or labels for your data
#
# e.g. hs1 10 20 0.5 paired=yes,special=no,myvar=0.5
#
# ordinarily, 'paired', 'special' and 'myvar' would not be allowed
restrict_parameter_names = no

# Unless set to 'yes', parameter names will be converted to lowercase
case_sensitive_parameter_names = no

# Calculate track statistics? Gives access to 
# min 
# max 
# average 
# n - number of data points in track
# sd - standard deviation
# var - variance (sd^2)
#
# Once these are calculated, you can refer to the value using plot(STATISTIC)
# e.g.
# <rule>
# condition = plot(average) < 0.5
# ...
calculate_track_statistics = yes

# The location of configuration and data files will be guessed if 
# (a) configuration file is not specified
# (b) data file paths are relative
# Circos will look in the following locations, where 
# SCRIPTPATH is the location of the 'circos' script (e.g. /usr/local/bin) and
# CWD is the current directory (where the 'circos' command was executed).
# All paths under CWD will be scanned first, then under SCRIPTPATH.
#
# [CWD,SCRIPTPATH]/.
# [CWD,SCRIPTPATH]/..
# [CWD,SCRIPTPATH]/etc/
# [CWD,SCRIPTPATH]/../etc
# [CWD,SCRIPTPATH]/../../etc
# [CWD,SCRIPTPATH]/data
# [CWD,SCRIPTPATH]/../data
# [CWD,SCRIPTPATH]/../../data
#
# If you would like to prepend this list with custom directories for
# data files, enter them as a CSV list here
# data_path = /home/martink/circos-tutorials 

# If the cache is static, it will always be used and will not be updated
# unless it is deleted (use -color_cache_rebuild on the command line).
# Otherwise, the cache will be updated if 
#  - config file is newer than cache file
#  - list of colors in config file is different than in cache file
color_cache_static = yes
color_cache_file   = circos.colorlist
color_lists_use    = yes
# if the directory is not defined, then the system will guess a temporary
# directory compatible with your operating system (using File::Temp)
# color_cache_dir    = /tmp

# Make some functions faster. This should always be 'yes' unless you
# want things to run slowly or suspect deep issues.
memoize = yes

# This is a debugging flag and should be set to 'no' for regular use
quit_on_dump = yes

offsets = 0,0

# Maximum number of image and data elements. If these are exceeded,
# Circos will quit with an error. These values are arbitrary, but in
# my experience images with significantly more data points than this
# are uninterpretable.

max_ticks            = 5000
max_ideograms        = 200
max_links            = 25000
max_points_per_track = 25000

# What to do when data is found for an ideogram that does not appear in the karyotype file.
# Set to 'skip' or 'exit'
undefined_ideogram = skip

# Number of iterations for determining ideogram sizes when
# relative scale is used.
relative_scale_iterations = 10

# min, max, average, mode             - based on scale statistics of ALL ideograms
# minadj, maxadj, averageadj, modeadj - based on scale statistics of adjacent ideograms
# 
# You can specify a fixed scale for spacing using a floating value
#
# e.g. relative_scale_spacing = 1.5
relative_scale_spacing    = mode

# What to do with out-of-range data, or data that extends beyond
# ideogram range. This applies to all plot types and links.  Set to
# 'hide' or 'trim' or 'fatal'. For 'hide' and 'trim', optionally add
# 'warn' to produce a warning. Data points that fall entirely
# outside of their ideogram will never be shown, regardless whether
# you've set 'hide' or 'trim'.
#
# This check DOES NOT trigger for data that falls on ideograms that
# are not being displayed.
#
# e.g.
# hide
# hide,warn
# trim
# trim,warn
# fatal
data_out_of_range = trim,warn # 

# Track default directory
track_defaults = etc/tracks

# Use round brushes for elements with thickness greater than round_brush_min_thickness?
round_brush_use           = yes
round_brush_min_thickness = 5

# Use anti aliasing, where possible? I've seen bugs in some gd libraries
# that cause artefacts to appear when lines are anti-aliased. If your
# image contains unexpected elements, turn aa off.
anti_aliasing = yes

# A parameter that must be set. Checks whether this file was imported.
housekeeping = yes

# Should parameter values be evaluated as code by default? This
# allows you to avoid eval() in expressions. However, this is
# dangerous because it does not return errors.
auto_eval = no
'''.format(font_color = font_color))
    file.close()

    file = open("conf/image.conf", "w")
    file.write('''<<include image.generic.conf>>
background = {}
'''.format(background_color))
    file.close()

    file = open("conf/ticks.conf", "w")
    file.write('''# TICKS
show_ticks          = yes
show_tick_labels    = yes
<ticks>
skip_first_label     = no
skip_last_label      = no
radius           = dims(ideogram,radius_outer)
thickness        = 2p
tick_separation      = 2p
min_label_distance_to_edge = 0p
label_separation = 5p
label_offset     = 5p
multiplier = 0.000001
color = {font_color}
<tick>
size     = 7p
thickness      = 1p
spacing        = 1000u
show_label     = no
</tick>
<tick>
size     = 10p
thickness      = 3p
spacing        = 5000u
show_label     = yes
suffix = " Mb"
label_size     = 20p
#format         = %s
</tick>
</ticks>
'''.format(font_color = font_color))
    file.close()
    return

if __name__ == "__main__":
    parser = ap.ArgumentParser()

    limit_group = parser.add_argument_group("limits")
    limit_group.add_argument("--content_min", "--min_GC_content", type=str, help="Minimum GC content", required=False, default = "0")
    limit_group.add_argument("--content_max", "--max_GC_content", type=str, help="Maximum GC content", required=False, default = "100")
    limit_group.add_argument("--skew_min", "--min_GC_skew", type=str, help="Minimum GC skew", required=False, default = "-1")
    limit_group.add_argument("--skew_max", "--max_GC_skew", type=str, help="Maximum GC skew", required=False, default = "1")

    color_group = parser.add_argument_group("colors")
    color_group.add_argument("-cc", "--GC_content_color", type=str, help="Color for GC content, in R, G, B format. Default: '23, 0, 115'", default = "23, 0, 115")
    color_group.add_argument("-sc", "--GC_skew_color", type=str, help="Color scheme for GC skew. For details on this, please read CIRCOS documentation. Default: 'eval(sprintf(\"rdbu-7-div-%%d\",remap_int(var(value),0,0,7,5)))'", default = 'eval(sprintf("rdbu-7-div-%d",remap_int(var(value),0,0,7,5)))')
    color_group.add_argument("-pc", "--CDS_positive_color", type=str, help="Color for positive CDSs, in R, G, B format. Default: '180, 205, 222'", default = '180, 205, 222')
    color_group.add_argument("-nc", "--CDS_negative_color", type=str, help="Color for negative CDSs, in R, G, B format. Default: '53, 176, 42'", default = '53, 176, 42')

    args = parser.parse_args()
    create_conf({"min_GC_content": args.min_GC_content,
                 "max_GC_content": args.max_GC_content,
                 "min_skew": args.min_GC_skew,
                 "max_skew": args.max_GC_skew},
                 args.GC_content_color,
                 args.GC_skew_color,
                 args.CDS_positive_color,
                 args.CDS_negative_color)
