
process STAR_INDEX {

    tag "${genome}_STAR_index"

    label 'high'

    publishDir "$params.star_index_path", mode : 'copy'
    
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
    fastq_command = params.single_end ? "--readFilesIn ${fastq}" : "--readFilesIn ${fastq[0]} ${fastq[1]}"
    readFilesCommand = params.readFilesCommand ? "--readFilesCommand zcat" : ""
    alignIntronMax_command = params.alignIntronMax ? "--alignIntronMax ${params.align_intron_max}" : "--alignIntronMax 10000"
    outFilterMismatchNoverReadLmax_command = params.outFilterMismatchNoverReadLmax ? "--outFilterMismatchNoverReadLmax ${outFilterMismatchNoverReadLmax}" : "--outFilterMismatchNoverReadLmax 0.04"
    outSAMtype_command = params.outSAMtype ? "--outSAMtype ${params.outSAMtype}" : "--outSAMtype BAM SortedByCoordinate"
    outFilterType_command = params.outFilterType ? "--outFilterType ${params.outFilterType}" : "--outFilterType BySJout"
    outFilterMultimapNmax_command = params.outFilterMultimapNmax ? "--outFilterMultimapNmax ${params.outFilterMultimapNmax}" : ""
    outFilterIntronMotifs_command = params.outFilterIntronMotifs ? "--outFilterIntronMotifs ${params.outFilterIntronMotifs}" : "--outFilterIntronMotifs RemoveNoncanonical"
    
    """
    #!/bin/bash
    
    source activate rnaseq

    STAR --genomeDir ${star_idx} \\
         --runThreadN ${task.cpus} \\
         ${outFilterMultimapNmax_command} \\
         ${outFilterType_command} \\
         ${gzip_command} \\
         ${outFilterMismatchNoverReadLmax_command} \\
         ${align_intron_max_command} \\
         ${fastq_command} \\
         ${outSAMtype_command} \\
         ${outFilterIntronMotifs_command} \\
         --outFileNamePrefix ${sampleID}_
 
    samtools index -@ 12 ${sampleID}_Aligned.sortedByCoord.out.bam
    """
}

