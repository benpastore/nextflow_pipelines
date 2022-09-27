process FASTQC {

    tag "${sampleID}_FASTQC"

    label 'medium'

    publishDir "$params.results/fastqc", mode : 'copy', pattern : '*_fastqc.{zip,html}'

    input : 
        tuple val(sampleID), val(fastq)
    
    output :
        path("*_fastqc.{zip,html}")

    script :
    fastq_command = params.single_end ? "${fastq}" : "${fastq[0]} ${fastq[1]}"
    """
    #!/bin/bash

    source activate rnaseq

    fastqc ${fastq_command} -o \$PWD

    """
}