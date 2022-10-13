process SALMON_INDEX {

    tag "${transcripts}_SALMON_index"

    label 'medium'

    publishDir "$params.salmon_index_path", mode : 'copy'

    input : 
        val transcripts
    
    output : 
        path("*")
        val("${params.salmon_index_path}"), emit : salmon_idx_ch
    
    script : 
    """
    #!/bin/bash

    source activate rnaseq

    salmon index --threads ${task.cpus} -t ${transcripts} --index \$PWD

    """

}

process SALMON {

    tag "${sampleID}_SALMON_map"

    label 'high'

    publishDir "$params.results/counts/salmon", mode : 'copy', pattern : "*.counts"

    input : 
        val salmon_idx
        tuple val(sampleID), val(fastq)

    output : 
        tuple val(sampleID), path("${sampleID}.salmon.counts"), emit : salmon_quant_ch
    
    script : 
    fastq_command = params.single_end ? "-r ${fastq}" : "-1 ${fastq[0]} -2 ${fastq[1]}"
    """
    #!/bin/bash

    source activate rnaseq

    salmon quant \\
        -i ${salmon_idx} \\
        --threads ${task.cpus} \\
        --libType SR \\
        ${fastq_command} \\
        -o ${sampleID}
    
    cat ${sampleID}/quant.sf | cut -f1,5 > ${sampleID}.salmon.counts
    """
}