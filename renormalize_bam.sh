#!/bin/bash

fdir=""

for f in $(ls $fdir/*.bam);
do

    #1. calculate normalization factor
    # print normalization factor to a file for reference in future 

    #2. bedtools genomecov

    #3. sort bedgraph

    #4. bedgraph --> bw

done