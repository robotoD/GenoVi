# GenoVi: Genome Visualizer Software
------
**GenoVi** is a visualization software tool for high quality circular representation of complete or draft genomes.
  
## Installation

GenoVi dependencies can be installed creating the following bioconda environment

```
conda create -n genovi python=3.7 deepnog circos 
```


 
## Dependencies
* [Circos](http://www.circos.ca/software/ "Circos")
* [DeepNOG](https://github.com/univieCUBE/deepnog "DeepNOG")
* Python 3.6.2 or later
* Perl

## Test tutorial

Type `GenoVi.py -h`

## Output files 
| Extension| Description|
| :-------------: |-------------:|
| .faa     |Protein FASTA file obtained from input .gbk file|
| bands.kar      | |
| .txt | input files for Circos visualization|
| prediction_deepnog.csv | COG prediction using the COG 2020 database by DeepNOG|


## Citation and License (CC BY-NC-SA? CC BY-NC? CC BY-NC-ND?)

GenoVi is under a BY-NC-SA Creative Commons License, Please cite.
Cumsille et al., 2021 
