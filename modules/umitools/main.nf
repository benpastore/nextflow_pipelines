process UMITOOLS_EXTRACT {

    label 'umitools'

    input : 
        tuple val(sampleID), val(fastqs)
    
    output : 
        tuple val(sampleID), path("*1.umi.fastq.gz"), path("*2.umi.fastq.gz"), emit : umi_extracted_fq

    script : 
    """
    #!/bin/bash

    R1=\$(basename ${fastqs[0]} .fastq.gz)
    R2=\$(basename ${fastqs[1]} .fastq.gz)
    
    umi_tools extract \\
        --bc-pattern=NNN \\
        --bc-pattern2=NNN \\
        --stdin=${fastqs[0]} \\
        --read2-in=${fastqs[1]} \\
        --stdout=\${R1}.umi.fastq.gz \\
        --read2-out=\${R2}.umi.fastq.gz 

    """

}

process UMITOOLS_DEDUP { 

    label 'umitools'

    input : 
        tuple val(sampleID), val(bam), val(bai)
    
    output : 
        tuple val(sampleID), path("*.dedup.bam"), path("*.dedup.bai"), emit : umi_dedup_bam
    
    script : 
    """
    #!/bin/bash

    name=\$(basename ${bam} .bam)
    umi_tools dedup \\
        -I ${bam} \\
        -s \$name.dedup.bam \\
        --umi_tag=UB \\
        --method=adjacency
    
    samtools index -@ 12 \$name.dedup.bam

    """

}

