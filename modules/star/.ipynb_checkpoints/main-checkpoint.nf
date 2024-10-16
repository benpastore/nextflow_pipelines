
process STAR_INDEX {

    tag "${genome}_STAR_index"

    label 'high'

    //publishDir "$params.star_index_path", mode : 'copy'
    
    input : 
        val genome
        val gtf
        val genome_dir
    
    output : 
        path("*")
        val("${params.star_index_path}"), emit : star_idx_ch

    script : 
    """
    #!/bin/bash

    source activate rnaseq
    
    STAR --runMode genomeGenerate \\
         --genomeFastaFiles ${genome} \\
         --sjdbGTFfile ${gtf} \\
         --runThreadN ${task.cpus} \\
         --genomeDir ${genome_dir}
    """
}

process STAR_ALIGN {

    tag "${sampleID}_STAR_align"

    label 'high'

    publishDir "$params.results/star_alignments/bam", mode : 'copy', pattern : '*.bam'
    publishDir "$params.results/star_alignments/bam", mode : 'copy', pattern : '*.bai'
    publishDir "$params.results/star_alignments/logs", mode : 'copy', pattern : '*Log.final.out'

    input : 
        val star_idx
        val gtf
        tuple val(sampleID), val(fastq)
    
    output :
        tuple val(sampleID), path("${sampleID}.Aligned.sortedByCoord.out.bam"), path("${sampleID}.Aligned.sortedByCoord.out.bam.bai"), emit : star_bams_ch
        path("*Log.final.out")

    script : 
    fastq_command = params.single_end ? "--readFilesIn ${fastq}" : "--readFilesIn ${fastq[0]} ${fastq[1]}"
    readFilesCommand = params.readFilesCommand ? "--readFilesCommand ${params.readFilesCommand}" : ""
    alignIntronMax_command = params.alignIntronMax ? "--alignIntronMax ${params.alignIntronMax}" : "--alignIntronMax 10000"
    outFilterMismatchNoverReadLmax_command = params.outFilterMismatchNoverReadLmax ? "--outFilterMismatchNoverReadLmax ${params.outFilterMismatchNoverReadLmax}" : "--outFilterMismatchNoverReadLmax 0.04"
    outSAMtype_command = params.outSAMtype ? "--outSAMtype ${params.outSAMtype}" : "--outSAMtype BAM SortedByCoordinate"
    outFilterType_command = params.outFilterType ? "--outFilterType ${params.outFilterType}" : "--outFilterType BySJout"
    outFilterIntronMotifs_command = params.outFilterIntronMotifs ? "--outFilterIntronMotifs ${params.outFilterIntronMotifs}" : "--outFilterIntronMotifs RemoveNoncanonical"
    outFilterMultimapNmax_command = params.outFilterMultimapNmax ? "--outFilterMultimapNmax ${params.outFilterMultimapNmax}" : "--outFilterMultimapNmax 20"
    outSAMprimaryFlag_command = params.outSAMprimaryFlag ? "--outSAMprimaryFlag ${params.outSAMprimaryFlag}" : "--outSAMprimaryFlag AllBestScore"
    outFilterScoreMinOverLread_command = params.outFilterScoreMinOverLread ? "--outFilterScoreMinOverLread ${params.outFilterScoreMinOverLread}" : "--outFilterScoreMinOverLread 0.66"
    outFilterMatchNminOverLread_command = params.outFilterMatchNminOverLread ? "--outFilterMatchNminOverLread ${params.outFilterMatchNminOverLread}" : "--outFilterMatchNminOverLread 0.66"
    outFilterMatchNmin_command = params.outFilterMatchNmin ? "--outFilterMatchNmin ${params.outFilterMatchNmin}" : "--outFilterMatchNmin 0"
    
    """
    #!/bin/bash
    
    source activate rnaseq

    STAR --genomeDir ${star_idx} \\
         --runThreadN ${task.cpus} \\
         ${fastq_command} \\
         ${readFilesCommand} \\
         ${alignIntronMax_command} \\
         ${outFilterMismatchNoverReadLmax_command} \\
         ${outSAMtype_command} \\
         ${outFilterType_command} \\
         ${outFilterMultimapNmax_command} \\
         ${outFilterIntronMotifs_command} \\
         ${outSAMprimaryFlag_command} \\
         ${outFilterScoreMinOverLread_command} \\
         ${outFilterMatchNminOverLread_command} \\
         ${outFilterMatchNmin_command} \\
         --outFileNamePrefix ${sampleID}.
 
    samtools index -@ 12 ${sampleID}.Aligned.sortedByCoord.out.bam
    """
}

