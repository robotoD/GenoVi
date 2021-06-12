#! /bin/bash

for file in *.fasta;
do
    echo $file
    for ((n=0;n<3;n++));
    do
      time python3 /d/Drive/GSoC/GC_analysis/GC_analysis/GC_analysis.py -i $file -o $file -w 5 -s 5
    done
done
      