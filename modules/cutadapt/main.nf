process CUTADAPT_QUALITY_TRIM {

    label 'cutadapt'

    input : 
        tuple val(sampleID), val(fastqR1), val(fastqR2)
    
    output : 
        tuple val(sampleID), path("*1*xq20*fastq.gz"), path("*2*xq20*fastq.gz"), emit : fastqs

    script : 
    """
    #!/bin/bash

    R1=\$(basename ${fastqR1} .fastq.gz)
    R2=\$(basename ${fastqR2} .fastq.gz)
    
    cutadapt \\
        -q 20,20 \\
        -o ${fastqR1} \\
        -p ${fastqR2} \\
        \${R1}.xq20.fastq.gz \${R2}.xq20.fastq.gz
    """

}

process CUTADAPT_TRIM_3PRIME_ADAPTER {

    label 'cutadapt'

    input : 
        tuple val(sampleID), val(fastqR1), val(fastqR2)
    
    output : 
        tuple val(sampleID), path("*1*x5p*fastq.gz"), path("*2*xq20*fastq.gz"), emit : fastqs

    script : 
    """
    #!/bin/bash

    R1=\$(basename ${fastqR1} .fastq.gz)
    R2=\$(basename ${fastqR2} .fastq.gz)
    
    cutadapt \
        -g AACAACGAGAAGATCGATGA \\                    # for R1
        -G TCATCGATCTTCTCGTTGTT \\                    # reverse complement for R2
        -o ${fastqR1} \\
        -p ${fastqR2} \\
        \${R1}.x5PrimeAdapter.fastq.gz \${R2}.x5PrimeAdapter.fastq.gz
    """

}

process CUTADAPT_TRIM_5PRIME_ADAPTER {

    label 'cutadapt'

    input : 
        tuple val(sampleID), val(fastqR1), val(fastqR2)
    
    output : 
        tuple val(sampleID), path("*1*x5p*fastq.gz"), path("*2*xq20*fastq.gz"), emit : fastqs

    script : 
    """
    #!/bin/bash

    R1=\$(basename ${fastqR1} .fastq.gz)
    R2=\$(basename ${fastqR2} .fastq.gz)
    
    cutadapt \\
        -a GGCGTCGCCATATTCTACTT \\
        -o ${fastqR1} \\
        -p ${fastqR2} \\
        \${R1}.x3primeAdapter.fastq.gz \${R2}.x3PrimeAdapter.fastq.gz
    """

}