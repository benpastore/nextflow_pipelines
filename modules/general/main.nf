process CHIP_SAMPLES_SHEET {

    label 'local'

    publishDir "$params.results/samples", mode : 'copy', pattern : '*'

    input : 
        val(samples)

    output : 
        path("controls.csv"), emit : control_ch
        path("fastq.csv"), emit : fq_ch
        path("replicates.csv"), emit : replicates_ch

    shell : 
    """
    #!/bin/bash

    source activate rnaseq

    which bwa

    python3 ${params.bin}/process_design_input.py -input ${samples}

    """
}

process MRNA_SAMPLES_SHEET {

    label 'local'

    publishDir "$params.results/samples", mode : 'copy', pattern : '*'

    input : 
        val(samples)

    output : 
        path("fastq.csv"), emit : fq_ch
        path("replicates.csv"), emit : replicates_ch

    shell : 
    """
    #!/bin/bash

    source activate rnaseq

    which bwa

    python3 ${params.bin}/process_design_input_mRNA.py -input ${samples}

    """
}

process MASTER_TABLE {

    label 'low'

    publishDir "$params.results/counts", mode : 'copy', pattern : "*.tsv"

    input :
        val project_name
        val counts
        
    output : 
        path("*.tsv"), emit : master_tables_ch
    
    script :
    gene_ids = params.gene_ids ? "-a ${params.gene_ids}" : ""
    """
    #!/bin/bash

    source activate rnaseq

    python3 ${params.bin}/make_master.py \\
        -f "${counts}" \\
        -o ${project_name} \\
        ${gene_ids}

    """
}

process RBIND_COUNTS {

    label 'low'

    publishDir "$params.results/counts", mode : 'copy', pattern : "*.tsv"

    input :
        tuple val(sampleID), path(feature_counts), path(salmon)

    output :
        path("*.tsv"), emit : counts_ch

    script : 
    """
    #!/bin/bash

    source activate rnaseq

    cat ${salmon} ${feature_counts} | grep -v "Name" | grep -v "Geneid" >> ${sampleID}.tsv
    """

}