# GenoVi: Genome Visualizer Software

**GenoVi** GenoVi generates circular genome representations for complete or draft bacterial genomes. GenoVi uses a pipeline that compiles several python scripts to generate all needed files for Circos software. Optionally, GenoVi uses DeepNOG to annotate COG classifications.

## Requirements
Circos 0.69-8
Python 3.7 or later
DeepNog 1.2.3
NumPy 1.20.2
Pandas 1.2.4
Biopython 1.79
CairoSVG 2.5.2
Perl 5
List::MoreUtils (Perl library)

## Installation

GenoVi dependencies can be installed creating the following bioconda environment

```
conda create -n genovi python=3.7 circos 
```

GenoVi then can be installed using pip

```
pip install genovi 
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

Please type `GenoVi.py -h` for complete usage instructions. Anyway, here are a few examples:

```
genovi -i input_test/Corynebacterium_alimapuense_VA37.gbk --status draft --color velvet --background black --font_color white --cogs_unclassified
```
![Corynebacterium alimapuense VA37](output_test/Corynebacterium_alimapuense_VA37-no_cogs.svg "Corynebacterium_alimapuense_VA37")
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


## Citation and License

GenoVi is under a BY-NC-SA Creative Commons License, Please cite.
Cumsille et al., 2021 

You may remix, tweak, and build upon this work even for commercial purposes, as long as you credit this work and license your new creations under the identical terms. 
