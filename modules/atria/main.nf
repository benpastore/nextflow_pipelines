process ATRIA {

    tag "${sampleID}_ATRIA"

    label 'atria'
    
    //publishDir "$params.results/trimmed_fastq", mode : 'copy', pattern : '*.fq.gz'
    publishDir "$params.results/atria/${sampleID}", mode : 'copy', pattern : '!*.fq.gz'

    input : 
        tuple val(sampleID), val(fastq)

    output :
        tuple val(sampleID), path("*.fq.gz"), emit : fq_ch
        path("*")


    script :
    minquality = params.minquality ? "-q ${params.minquality}" : "-q 20"
    minlength = params.minlength ? "--length-range ${params.lenrange}" : "--length-range 50:500"
    fastq_command = params.single_end ? "-r ${fastq}" : "--read1 ${fastq[0]} --read2 ${fastq[1]}"
    """
    #!/bin/bash
    #singularity run docker://benpasto/atria:latest /run_atria -h
    /bin/run_atria ${minquality} ${minlength} -t ${task.cpus} ${fastq_command}

    """

}
