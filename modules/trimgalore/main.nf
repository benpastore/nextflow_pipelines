process TRIM_GALORE {

    tag "${sampleID}_TrimGalore"

    label 'medium'
    
    //publishDir "$params.results/trimmed_fastq", mode : 'copy', pattern : '*.fq.gz'
    publishDir "$params.results/trimmed_fastq", mode : 'copy', pattern : '*report.txt'

    input : 
        tuple val(sampleID), val(fastq)

    output :
        tuple val(sampleID), path("${sampleID}*.fq.gz"), emit : fq_ch

    script :
    trimn = params.trimN ? "--trim-n" : ""
    minquality = params.minquality ? "-q ${params.minquality}" : "-q 20"
    minlength = params.minlength ? "--length ${params.minlength}" : "--length 20"
    fastq_command = params.single_end ? "${fastq}" : "${fastq[0]} ${fastq[1]}"
    paired = params.single_end ? "" : "--paired"
    """
    #!/bin/bash
    
    source activate rnaseq

    trim_galore ${minquality} ${minlength} ${trimn} -j ${task.cpus} ${paired} --gzip --basename ${sampleID} ${fastq_command}

    """

}