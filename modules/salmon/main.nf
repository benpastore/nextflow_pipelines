process SALMON_INDEX {

    tag "${transcripts}_SALMON_index"

    label 'medium'

    //publishDir "$params.salmon_index_path", mode : 'copy'

    input : 
        val transcripts
        val genome_dir
    
    output : 
        path("${params.salmon_index_path}"), emit : salmon_idx_ch
        val genome_dir
    
    script : 
    additional_salmon_commands = params.salmon_commands ? "${params.salmon_commands}" : ""
    """
    #!/bin/bash

    source activate rnaseq

    salmon index --threads ${task.cpus} -t ${transcripts} ${additional_salmon_commands} --index ${params.salmon_index_path}

    [ ! -d ${genome_dir} ] && mkdir -p ${genome_dir}

    cp * ${genome_dir}

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
    salmon_libtype = params.salmon_libtype ? "${params.salmon_libtype}" : "ISR"
    """
    #!/bin/bash

    source activate rnaseq

    salmon quant \\
        -i ${salmon_idx} \\
        --threads ${task.cpus} \\
        --libType ${salmon_libtype} \\
        ${fastq_command} \\
        -o ${sampleID}
    
    cat ${sampleID}/quant.sf | cut -f1,5 > ${sampleID}.salmon.counts
    """
}