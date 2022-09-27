process MASTER_TABLE {

    label 'low'

    publishDir "$params.results/counts", mode : 'copy', pattern : "*.tsv"

    input :
        val project_name
        val counts
        val gene_ids
        
    output : 
        path("*.tsv"), emit : master_tables_ch
    
    script :
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