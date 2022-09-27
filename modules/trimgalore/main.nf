process TRIM_GALORE {

    tag "${sampleID}_TrimGalore"

    label 'medium'
    
    publishDir "$params.results/trimmed_fastq", mode : 'copy', pattern : '*.fq.gz'
    publishDir "$params.results/trimmed_fastq", mode : 'copy', pattern : '_trimming_report.txt'

    input : 
        tuple val(sampleID), val(fastq)

    output :
        tuple val(sampleID), path("${sampleID}*.fq.gz"), emit : fq_ch

    script :
    fastq_command = params.single_end ? "${fastq}" : "${fastq[0]} ${fastq[1]}"
    paired = params.single_end ? "" : "--paired"
    """
    #!/bin/bash
    
    source activate rnaseq

    trim_galore -j ${task.cpus} ${paired} --gzip --basename ${sampleID} ${fastq_command}

    """

}