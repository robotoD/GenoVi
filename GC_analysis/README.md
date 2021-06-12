[![Build status](https://travis-ci.org/tonyyzy/GC_analysis.svg?branch=master)](https://travis-ci.org/tonyyzy/GC_analysis)
[docs](https://gc-analysis.readthedocs.io/en/latest/index.html)
# GC-analysis
A command-line utility for calculating GC percentages of genome sequences

# Quick starter
Calculate the GC content of chromosome 17 of the human reference genome with window size (or span) = 5 and shift (or step) = 5. Input fasta file is `GRCh38-Chrom17.fasta` and output wiggle file is `GRCh38-Chrom17.wig`. Note that the output file's extension is added by the program.
```
~ $ GC_analysis -i GRCh38-Chrom17.fasta -w 5 -s 5 -o GRCh38-Chrom17
```

# Installation guide
Note that pyBigWig can only be used under linux environment. To work with Windows system, the Docker image can be used as shown below. Alternatively, you can clone the repository, comment out `import pyBigWig` and the script would work but without BigWig support.

1. Pip install GC_analysis (NB. Python3 is recommanded but GC_analysis should work with Python2 as well)
```
pip3 install GC_analysis
```
Then `GC_analysis.py` command will be available globally.
```
GC_analysis.py -i [INPUT] -o [OUTPUT] -w [window size] -s [shift]
```

2. Run the python script directly. Please ensure you have python3 installed with pyBigwig and Biopython.
Clone the github repository and install packages.
```
git clone https://github.com/tonyyzy/GC_analysis
cd GC_analysis
pip3 install -r requirements.txt
```
run the script from `GC_analysis` directory.
```
python3 ./GC_analysis/GC_analysis.py -i [INPUT] -o [OUTPUT] -w [window size] -s [shift]
```

3. Use the packaged binary.
```
mkdir ~/GC_analysis
cd ~/GC_analysis
wget https://github.com/tonyyzy/GC_analysis/releases/download/v0.3/GC_analysis
```
Execute the binary command
```
GC_analysis -i [INPUT] -o [OUTPUT] -w [window size] -s [shift]
```

4. Use the Docker image.
Firstly, pull the docker image (around 384 MB)
```
docker pull tonyyzy/gc_analysis
```
To use input files outside the container and save output files on your computer, the `-v` volume mapping option will be used. You will need to know the absolute path of the directory you want to map (which can be found out with `pwd`).
```
docker run -v /your/local/path:/app tonyyzy/gc_analysis GC_analysis -i /app/yours.fasta -o /app/yours -w 5 -s 5
```
This option maps `/your/local/path` to `/app` under the container's root directory. Your result file will be saved to `/your/local/path/yours.wig`.

# Command-line options
```
~ $ GC_analysis -h
usage: GC_analysis [-h] -i INPUT_FILE -w WINDOW_SIZE -s SHIFT [-o OUTPUT_FILE]
                   [-ot] [-f {wiggle,gzip,bigwig}]

required named arguments:

-i INPUT_FILE, --input_file INPUT_FILE
INPUTFILE: Name of the input file in FASTA format

-w WINDOW_SIZE, --window_size WINDOW_SIZE
WINDOW_SIZE: Number of base pairs that the GC percentage is calculated for

-s SHIFT, --shift SHIFT
SHIFT: The shift increment (step size)

optional arguments:

-h, --help
Show the help message and exit

-o OUTPUT_FILE, --output_file OUTPUT_FILE
OUTPUT_FILE: Name of the output file

-ot, --omit_tail
Use if the trailing sequence should be omitted. Default behaviour is to retain the leftover sequence.

-f {wiggle,bigwig,gzip}, --output_format {wiggle,bigwig,gzip}
Choose output formats from wiggle, bigwig or gzip compressed wiggle file.

-one, --one_file
Force one file output

```
## Example usage
1. Calculate the GC content of chromosome 17 of the human reference genome, the percentage is calculated over five base pairs (window_size), and the window is shifted by five base pairs every time (i.e. there is no overlapping base paires in each entry).
```
~ $ GC_analysis -i GRCh38-Chrom17.fasta -w 5 -s 5 -o GRCh38-Chrom17
```

2. By default, the GC percentage of the trailing sequence is calculated and appended to the end of the output file. For example, with the following input
```
~ $ GC_analysis -i examaple1.fasta -w 5 -s 5 -o with_tail
```
and `example1.fasta` is
```
>chr1
AAAAACC
```
the generated `with_tail.wig` will look like
```
track type=wiggle_0 name="GC percentage" description="chr1"
variableStep chrom=chr1 span=5
1	0
6	100
```
If it is desirable to omit the trailing sequence in the result, the `-ot` or `--omit_tail` option can be used. For example
```
~ $ GC_analysis -i examaple1.fasta -w 5 -s 5 -o without_tail -ot
```
will generate output file `without_tail` with the following content
```
track type=wiggle_0 name="GC percentage" description="chr1"
variableStep chrom=chr1 span=5
1	0
```

3. The program support three output file formats, wiggle, bigwig and gzip compressed wiggle file.
Wiggle output file follows the [UCSC variableStep format definition](https://genome.ucsc.edu/goldenpath/help/wiggle.html). Wiggle file is the default output format. The output format can be changed with `-f` or `--format` option.
```
~ $ GC_analysis -i GRCh38-Chrom17.fasta -w 5 -s 5 -o GRCh38-Chrom17
```
and
```
~ $ GC_analysis -i GRCh38-Chrom17.fasta -w 5 -s 5 -o GRCh38-Chrom17 -f wiggle
```
will generate `GRCh38-Chrom17.wig` as the output file.

```
~ $ GC_analysis -i GRCh38-Chrom17.fasta -w 5 -s 5 -o GRCh38-Chrom17 -f gzip
```
will generate `GRCh38-Chrom17.wig.gz` as the output file. Decompress `GRCh38-Chrom17.wig.gz` will give you the same wiggle file as choosing wiggle as the output format.

```
~ $ GC_analysis -i GRCh38-Chrom17.fasta -w 5 -s 5 -o GRCh38-Chrom17 -f bigwig
```
will generate `GRCh38-Chrom17.bw` as the output file. It should be noted that bigwig format does not allow overlapping bases, which means that `-w 5 -s 3` is an invalid option with choosing bigwig as the output format. In this case, where shift is smaller than window size and bigwig format is specified, the program will generate a wiggle file instead and output a warning message.

```
~ $ GC_analysis -i GRCh38-Chrom17.fasta -w 5 -s 3 -o GRCh38-Chrom17 -f bigwig
WARNING! BigWig file does not allow overlapped items. A wiggle file was generated instead.
```

4. If an output filename is not given, the result will be written to stdout. If the output filename is not given and a file format other than wiggle was chosen, the program will automatically output the result to stdout and give you a warning before and after the result.
Eg. 
```
GC_analysis -i example1.fasta -w 5 -s 3 -f bigwig
WARNING! BigWig file does not allow overlapped items. A wiggle file will be generated instead.
WARNING! An output filename is needed to save output as bigwig. The result is shown below:
track type=wiggle_0 name="GC percentage" description="chr1"
variableStep chrom=chr1 span=5
1       0
4       50
WARNING! BigWig file does not allow overlapped items. A wiggle file was generated instead.
WARNING! An output filename is needed to save output as bigwig. The result is shown above.
```
5. If the input FASTA file contains multiple sequences and the results should be written to a single file instead of one sequence per file, you can use the `-one` optional arguments.
```
GC_analysis -i multiple.fasta -o multiple -w 5 -s 5 -one
```
If `-one` is not specified, each sequence's GC result will be written to one file. The filenames will be the given filename followed by \"_seq\" + the sequence's number.
## Timing againts human chromosomes
<details><summary><b>Click for raw data table</b></summary>
<p>
<table class="tg">
  <tr>
    <th class="tg-us36">Entry</th>
    <th class="tg-us36">Human chromosome</th>
    <th class="tg-us36">No. of base pairs</th>
    <th class="tg-us36">Average real time - single thread (s)</th>
    <th class="tg-us36">Average real time - multi threads (s)</th>
  </tr>
  <tr>
    <td class="tg-us36">CM000663.2.fasta</td>
    <td class="tg-us36">1</td>
    <td class="tg-us36">248956422</td>
    <td class="tg-us36">288.429</td>
    <td class="tg-us36">179.221</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000664.2.fasta</td>
    <td class="tg-us36">2</td>
    <td class="tg-us36">242193529</td>
    <td class="tg-us36">276.355</td>
    <td class="tg-us36">169.611</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000665.2.fasta</td>
    <td class="tg-us36">3</td>
    <td class="tg-us36">198295559</td>
    <td class="tg-us36">227.528</td>
    <td class="tg-us36">135.637</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000666.2.fasta</td>
    <td class="tg-us36">4</td>
    <td class="tg-us36">190214555</td>
    <td class="tg-us36">217.846</td>
    <td class="tg-us36">153.091</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000667.2.fasta</td>
    <td class="tg-us36">5</td>
    <td class="tg-us36">181538259</td>
    <td class="tg-us36">205.623</td>
    <td class="tg-us36">123.858</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000668.2.fasta</td>
    <td class="tg-us36">6</td>
    <td class="tg-us36">170805979</td>
    <td class="tg-us36">193.209</td>
    <td class="tg-us36">117.180</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000669.2.fasta</td>
    <td class="tg-us36">7</td>
    <td class="tg-us36">159345973</td>
    <td class="tg-us36">183.445</td>
    <td class="tg-us36">109.135</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000670.2.fasta</td>
    <td class="tg-us36">8</td>
    <td class="tg-us36">145138636</td>
    <td class="tg-us36">166.607</td>
    <td class="tg-us36">98.632</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000671.2.fasta</td>
    <td class="tg-us36">9</td>
    <td class="tg-us36">138394717</td>
    <td class="tg-us36">157.142</td>
    <td class="tg-us36">93.898</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000672.2.fasta</td>
    <td class="tg-us36">10</td>
    <td class="tg-us36">133797422</td>
    <td class="tg-us36">150.872</td>
    <td class="tg-us36">92.371</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000673.2.fasta</td>
    <td class="tg-us36">11</td>
    <td class="tg-us36">135086622</td>
    <td class="tg-us36">154.003</td>
    <td class="tg-us36">92.498</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000674.2.fasta</td>
    <td class="tg-us36">12</td>
    <td class="tg-us36">133275309</td>
    <td class="tg-us36">150.533</td>
    <td class="tg-us36">90.807</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000675.2.fasta</td>
    <td class="tg-us36">13</td>
    <td class="tg-us36">114364328</td>
    <td class="tg-us36">129.951</td>
    <td class="tg-us36">77.498</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000676.2.fasta</td>
    <td class="tg-us36">14</td>
    <td class="tg-us36">107043718</td>
    <td class="tg-us36">121.008</td>
    <td class="tg-us36">71.970</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000677.2.fasta</td>
    <td class="tg-us36">15</td>
    <td class="tg-us36">101991189</td>
    <td class="tg-us36">115.194</td>
    <td class="tg-us36">68.336</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000678.2.fasta</td>
    <td class="tg-us36">16</td>
    <td class="tg-us36">90338345</td>
    <td class="tg-us36">103.169</td>
    <td class="tg-us36">60.799</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000679.2.fasta</td>
    <td class="tg-us36">17</td>
    <td class="tg-us36">83257441</td>
    <td class="tg-us36">94.353</td>
    <td class="tg-us36">55.729</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000680.2.fasta</td>
    <td class="tg-us36">18</td>
    <td class="tg-us36">80373285</td>
    <td class="tg-us36">92.020</td>
    <td class="tg-us36">53.395</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000681.2.fasta</td>
    <td class="tg-us36">19</td>
    <td class="tg-us36">58617616</td>
    <td class="tg-us36">67.506</td>
    <td class="tg-us36">39.308</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000682.2.fasta</td>
    <td class="tg-us36">20</td>
    <td class="tg-us36">64444167</td>
    <td class="tg-us36">74.048</td>
    <td class="tg-us36">43.280</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000683.2.fasta</td>
    <td class="tg-us36">21</td>
    <td class="tg-us36">46709983</td>
    <td class="tg-us36">53.633</td>
    <td class="tg-us36">31.118</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000684.2.fasta</td>
    <td class="tg-us36">22</td>
    <td class="tg-us36">50818468</td>
    <td class="tg-us36">57.466</td>
    <td class="tg-us36">33.701</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000685.2.fasta</td>
    <td class="tg-us36">X</td>
    <td class="tg-us36">156040895</td>
    <td class="tg-us36">176.895</td>
    <td class="tg-us36">105.408</td>
  </tr>
  <tr>
    <td class="tg-us36">CM000686.2.fasta</td>
    <td class="tg-us36">Y</td>
    <td class="tg-us36">57227415</td>
    <td class="tg-us36">67.016</td>
    <td class="tg-us36">38.142</td>
  </tr>
  <tr>
    <td class="tg-us36">J01415.2.fasta</td>
    <td class="tg-us36">MT</td>
    <td class="tg-us36">16569</td>
    <td class="tg-us36">0.231</td>
    <td class="tg-us36">0.397</td>
  </tr>
</table>
  
  </p>
</details>




![Execution time vs. number of base pairs plot](https://github.com/tonyyzy/GC_analysis/blob/master/tests/time_profile/GC_time_profile.png "execution time plot")
\* 1) Real time data is the average of three runs; 2) GC_analysis parameters for each run is `-w 5 -s 5`; 3) `Serial` data is collected with the `Master` branch, `Parallel` data is collected with the `Parallel` branch.

As can be seen from the plot, `GC_analysis` scales well with number of base pairs, resulted a linear relationship between the execution time and the size of the chromosomes. Although multi-threaded version can provide ~1.7x speed improvement, it has a significantly higher memory consumption, hence it's not recommended.


## (EXPERIMENTAL!!!) Multi-threaded GC_analysis
Git clone the `parallel` branch from GitHub repo:
```
git clone --single-branch -b parallel https://github.com/tonyyzy/GC_analysis
```
Execute as normal from `GC_analysis` directory
```
~ python3 ./scripts/GC_analysis.py -i GRCh38-Chrom17.fasta -w 5 -s 5 -o GRCh38-Chrom17
```

This multithreading implementation is a very crude one and only result in ~1.7x speed up. A large amount of RAM is needed to store out-of-order intermediate results for sorting.
