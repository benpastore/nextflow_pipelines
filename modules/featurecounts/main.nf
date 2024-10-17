process FEATURECOUNTS {

    tag "${sampleID}_FeatureCounts"

    label 'medium'

    publishDir "$params.results/counts/feature_counts", mode: 'copy', pattern : "*.txt"
    publishDir "$params.results/counts/feature_counts", mode: 'copy', pattern : "*.summary"
    publishDir "$params.results/counts/feature_counts", mode: 'copy', pattern : "*.log"
    publishDir "$params.results/counts/feature_counts", mode: 'copy', pattern : "*.tsv"

    input :
        tuple val(sampleID), path(bam), path(bai)
        val gtf

    output :
        tuple val(sampleID), path("${sampleID}.counts.tsv"), emit : feature_counts_ch
        path("${sampleID}.feature_counts.txt")
        path("*.feature_counts.txt.summary")
        path("*.feature_counts.log")
    
    script :
    libtype_command = params.single_end ? "" : "-p"
    orientation_command = params.reverse_stranded ? "-s 2" : "-s 0"
    gene_id_command = params.identifier ? "-g ${params.identifier}" : "-g gene_id"
    feature_command = params.feature ? "-t ${params.feature}" : "-t exon"
    multimap_command = params.count_multimappers ? "-M" : ""
    fraction_command = params.fraction_counts ? "--fraction" : ""
    multioverlap_command = params.multi_overlap ? "-O" : ""
    """
    #!/bin/bash

    source activate rnaseq

    featureCounts \\
        ${feature_command} \\
        ${gene_id_command} \\
        ${orientation_command} \\
        ${libtype_command} \\
        ${multimap_command} \\
        ${fraction_command} \\
        ${multioverlap_command} \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -o ${sampleID}.feature_counts.txt \\
        ${bam} 2> ${sampleID}.feature_counts.log

    cat ${sampleID}.feature_counts.txt | grep -v "^#" | cut -f 1,7 | sed -e 's/.umi.dedup.bam//g' | sed -e 's/_Aligned.sortedByCoord.out.bam//g' | awk '(NR>1)' > ${sampleID}.counts.tsv

    cp ${sampleID}.feature_counts.txt.summary tmp.summary 
    rm ${sampleID}.feature_counts.txt.summary
    cat tmp.summary | sed -e 's/.umi.dedup.bam//g' | sed -e 's/_Aligned.sortedByCoord.out.bam//g' > ${sampleID}.feature_counts.txt.summary

    """

}