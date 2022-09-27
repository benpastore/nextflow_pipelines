
process STAR_INDEX {

    tag "${genome}_STAR_index"

    label 'high'

    publishDir "$params.star_index", mode : 'copy'
    
    input : 
        val genome
        val gtf
    
    output : 
        path("*")
        val("${params.star_index}"), emit : star_idx_ch

    script : 
    """
    #!/bin/bash

    source activate rnaseq
    
    STAR --runMode genomeGenerate \\
         --genomeDir \$PWD \\
         --genomeFastaFiles ${genome} \\
         --sjdbGTFfile ${gtf} \\
         --runThreadN ${task.cpus}
    """
}

process STAR_ALIGN {

    tag "${sampleID}_STAR_align"

    label 'high'

    publishDir "$params.results/star_alignments/bam", mode : 'copy', pattern : '*.bam'
    publishDir "$params.results/star_alignments/bam", mode : 'copy', pattern : '*.bai'
    publishDir "$params.results/star_alignments/logs", mode : 'copy', pattern : '*_Log.final.out'

    input : 
        val star_idx
        val gtf
        tuple val(sampleID), val(fastq)
    
    output :
        tuple val(sampleID), path("${sampleID}_Aligned.sortedByCoord.out.bam"), path("${sampleID}_Aligned.sortedByCoord.out.bam.bai"), emit : star_bams_ch
        path("*_Log.final.out")

    script : 
    fastq_command = params.single_end ? "${fastq}" : "${fastq[0]} ${fastq[1]}"
    """
    #!/bin/bash
    
    source activate rnaseq

    STAR --genomeDir ${star_idx} \\
         --runThreadN ${task.cpus} \\
         --outFilterMultimapNmax 999 \\
         --outFilterType BySJout \\
         --readFilesCommand zcat \\
         --alignIntronMax 10000 \\
         --outFileNamePrefix ${sampleID}_ \\
         --outFilterMismatchNoverReadLmax 0.04 \\
         --readFilesIn ${fastq_command} \\
         --outSAMtype BAM SortedByCoordinate \\
         --outFilterIntronMotifs RemoveNoncanonical
    
    samtools index -@ 12 ${sampleID}_Aligned.sortedByCoord.out.bam
    """
}

