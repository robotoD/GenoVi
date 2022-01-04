# GenoVi: Genome Visualizer Software

**GenoVi** generates circular genome representations for complete or draft bacterial and archaeal genomes. **GenoVi** pipeline combines several python scripts to automatically generate all needed files for Circos, including customizable options for color palettes, fonts, font format, background color and scaling options for complete genomes comprising more than 1 replicon. Optionally, **GenoVi** built-in workflow integrates DeepNOG to annotate COG categories using alignment-free methods with user-defined thresholds.

![Diagram](Figures/Diagram.png "Diagram")

## Installation

GenoVi dependencies can be installed creating the following bioconda environment

```
conda create -n genovi python=3.7 circos 
```
Activate the environment
```
conda activate genovi
```
GenoVi can then be installed using pip

```
pip install genovi=0.2.0
```
 
## Dependencies
* [Circos](http://www.circos.ca/software/ "Circos")
* [DeepNOG](https://github.com/univieCUBE/deepnog "DeepNOG")
* Python 3.7 or later
* Perl
* NumPy 1.20.2
* Pandas 1.2.4 
* Biopython 1.79
* CairoSVG 2.5.2 

## Usage

Please type `genovi -h` for complete usage. 

```
genovi [-h] [options ..] -i input_file -s status
```

* -i, --input_file. GenBank input file path.
* -o, --output_file.  Output file name. Default: circos.
* -s, --status. Draft or complete. Complete genomes are drawn as separate circles for each contig.

### Information:
* -h, --help. Shows this help message and exit
* --version. Shows the currently installed version of genovi.

### COGs:
* -cu, --cogs_unclassified. Do not classify each coding sequence into Clusters of Orthologous Groups of proteins (COGs)
* -b, --deepnog_confidence_threshold. DeepNOG confidence threshold range [0,1] Default: 0. If provided, predictions below the threshold are discarded.

### Format:
* -a, --alignment. When a --status complete is specified, this flag defines the alignment of each individual contig. Options: center, top, bottom, A (First on top), < (first to the left), U (Two on top, the rest below). By default this is defined by contig sizes
* --scale. When using --status complete, whether to use a different scale format to ensure visibility. Options: variable, linear, sqrt. Default: sqrt
* -k, --keep_temporary_files. Keep temporary files
* -w, --window. Window size (base pair) to assign a GC analysis. Default: 5000
* -v, --verbose. Verbose

### Text:
* -c, --captions_not_included.  Do not include captions in the figure.
* -cp, --captions_position. Captions position, options: left, right, auto.
* -t, --title. Figure 's title.
* --title_position. Title position {center,top,bottom}
* --italic_words. How many words should be written in italic in the title. Default: 2
* --size. Display genome size of each independent circular representation.

### Colors:
* -cs, --color_scheme. Prebuilt color scheme to use for CDS, RNAs and GC analysis. Options: neutral, blue, purple, soil, grayscale, velvet, pastel, ocean, wood, beach, desert, ice, island, forest, toxic, fire, spring.
* -bc, --background. Background color. Default: transparent
* -fc, --font_color. Font color. Default: black
* -pc, --CDS_positive_color. Color for positive CDSs, in R, G, B format. Default: '180, 205, 222'
* -nc, --CDS_negative_color. Color for negative CDSs, in R, G, B format. Default: '53, 176, 42'
* -tc, --tRNA_color. Color for tRNAs, in R, G, B format. Default: '150, 5, 50'
* -rc, --rRNA_color. Color for rRNAs, in R, G, B format. Default: '150, 150, 50'
* -cc, --GC_content_color. Color for GC content, in R, G, B format. Default: '23, 0, 115'
* -sc, --GC_skew_color- Color scheme for positive and negative GC skew. A pair of RGB colors. Default: '140, 150, 198 - 158, 188, 218'
* -sl, --GC_skew_line_color. Color for GC skew line. Default: black

### Prebuilt color palettes
![Prebuilt color palettes](Figures/Palettes.png "Prebuilt color palettes")


## Tutorial

### Draft circular genome representation without COGs classification
```
genovi -i input_test/Corynebacterium_alimapuense_VA37.gbk --status draft --color velvet --background "255,247,248" --cogs_unclassified
```
![Corynebacterium alimapuense VA37](output_test/Corynebacterium_alimapuense_VA37-no_cogs.svg "Corynebacterium_alimapuense_VA37")

### Complete circular genome representation
```
genovi -i input_test/Acinetobacter_radioresistens_DD78.gbk --status complete --color pastel
```
![Acinetobacter radioresistens DD78](output_test/Acinetobacter_radioresistens_DD78.svg "Acinetobacter_radioresistens_DD78")


## Output files 
| Extension| Description|
| :-------------: |-------------:|
| .svg     |Scalable vector graphic representing the genome.|
| .png     |Pixel-defined image representing the genome.|
| .faa     |Protein FASTA file obtained from input .gbk file|
| bands.kar      | Circos input file.|
| .txt | input files for Circos visualization|
| prediction_deepnog.csv | COG prediction using the COG 2020 database by DeepNOG.|

## Additional information
For further information, please read the user guide. 

## Citation and License

GenoVi is under a BY-NC-SA Creative Commons License, Please cite.
Cumsille et al., 2022 

You may remix, tweak, and build upon this work even for commercial purposes, as long as you credit this work and license your new creations under the identical terms. 
