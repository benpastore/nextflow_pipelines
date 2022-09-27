process COUNTER {

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
    libtype = params.single_end ? "" : "-p"
    """
    #!/bin/bash

    source activate rnaseq

    featureCounts \\
        -t exon \\
        -g gene_id \\
        -s 2 \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -o ${sampleID}.feature_counts.txt \\
        ${libtype} \\
        ${bam} 2> ${sampleID}.feature_counts.log

    cat ${sampleID}.feature_counts.txt | grep -v "^#" | cut -f 1,7 | sed -e 's/.umi.dedup.bam//g' | sed -e 's/_Aligned.sortedByCoord.out.bam//g' | awk '(NR>1)' > ${sampleID}.counts.tsv

    cp ${sampleID}.feature_counts.txt.summary tmp.summary 
    rm ${sampleID}.feature_counts.txt.summary
    cat tmp.summary | sed -e 's/.umi.dedup.bam//g' | sed -e 's/_Aligned.sortedByCoord.out.bam//g' > ${sampleID}.feature_counts.txt.summary

    """

}