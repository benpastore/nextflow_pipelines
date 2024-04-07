#!/usr/bin/env nextflow


test = Channel.of("./a.txt")

test
    .map{ it -> [it, it] }
    .map{ it -> $it[0]"\t"$it[1]}
    .collectFile( name : "test.tsv", storeDir : ".", newLine : true)