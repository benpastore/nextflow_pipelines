
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
        val("${genome_dir}"), emit : star_idx_ch

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
    publishDir "$params.results/star_alignments/counts", mode : 'copy', pattern : '*ReadsPerGene.out.tab'

    input : 
        val star_idx
        val gtf
        tuple val(sampleID), val(fastq)
    
    output :
        tuple val(sampleID), path("${sampleID}.Aligned.sortedByCoord.out.bam"), path("${sampleID}.Aligned.sortedByCoord.out.bam.bai"), emit : star_bams_ch // path("${sampleID}.Aligned.sortedByCoord.out.bam.bai"), emit : star_bams_ch
        path("*Log.final.out")
        tuple val(sampleID), path("${sampleID}.Aligned.sortedByCoord.out.bam"), path("${sampleID}.Aligned.sortedByCoord.out.bam.bai"), emit : star_bam_bai_ch
        path("*ReadsPerGene.out.tab")

    script : 
    fastq_command = params.single_end ? "--readFilesIn ${fastq}" : "--readFilesIn ${fastq[0]} ${fastq[1]}"
    readFilesCommand = params.readFilesCommand ? "--readFilesCommand ${params.readFilesCommand}" : ""
    alignIntronMax_command = params.alignIntronMax ? "--alignIntronMax ${params.alignIntronMax}" : "--alignIntronMax 10000"
    outFilterMismatchNoverReadLmax_command = params.outFilterMismatchNoverReadLmax ? "--outFilterMismatchNoverReadLmax ${params.outFilterMismatchNoverReadLmax}" : "--outFilterMismatchNoverReadLmax 0.04"
    outSAMtype_command = params.outSAMtype ? "--outSAMtype ${params.outSAMtype}" : "--outSAMtype BAM SortedByCoordinate"
    outFilterType_command = params.outFilterType ? "--outFilterType ${params.outFilterType}" : "--outFilterType BySJout"
    outFilterIntronMotifs_command = params.outFilterIntronMotifs ? "--outFilterIntronMotifs ${params.outFilterIntronMotifs}" : "--outFilterIntronMotifs RemoveNoncanonical"
    outFilterMultimapNmax_command = params.outFilterMultimapNmax ? "--outFilterMultimapNmax ${params.outFilterMultimapNmax}" : "--outFilterMultimapNmax 20"
    outSAMprimaryFlag_command = params.outSAMprimaryFlag ? "--outSAMprimaryFlag ${params.outSAMprimaryFlag}" : "--outSAMprimaryFlag OneBestScore"
    outFilterScoreMinOverLread_command = params.outFilterScoreMinOverLread ? "--outFilterScoreMinOverLread ${params.outFilterScoreMinOverLread}" : "--outFilterScoreMinOverLread 0.66"
    outFilterMatchNminOverLread_command = params.outFilterMatchNminOverLread ? "--outFilterMatchNminOverLread ${params.outFilterMatchNminOverLread}" : "--outFilterMatchNminOverLread 0.66"
    outFilterMatchNmin_command = params.outFilterMatchNmin ? "--outFilterMatchNmin ${params.outFilterMatchNmin}" : "--outFilterMatchNmin 0"
    outFilterMismatchNmax_command = params.outFilterMismatchNmax ? "--outFilterMismatchNmax ${params.outFilterMismatchNmax}" : "--outFilterMismatchNmax 10"
    outFilterScoreMin_command = params.outFilterScoreMin ? "--outFilterScoreMin ${params.outFilterScoreMin}" : "--outFilterScoreMin 0"
    outFilterMultimapScoreRange_command = params.outFilterMultimapScoreRange ? "--outFilterMultimapScoreRange ${params.outFilterMultimapScoreRange}" : "--outFilterMultimapScoreRange 1"
    alignEndsType_command = params.alignEndsType ? "--alignEndsType ${params.alignEndsType}" : "--alignEndsType Local"
    outSAMunmapped_command = params.outSAMunmapped ? "--outSAMunmapped ${params.outSAMunmapped}" : "--outSAMunmapped None"
    genomeLoad_command = params.genomeLoad ? "--genomeLoad ${params.genomeLoad}" : "--genomeLoad NoSharedMemory" 
    outReadsUnmapped_command = params.outReadsUnmapped ? "--outReadsUnmapped ${params.outReadsUnmapped}" :  "--outReadsUnmapped None" 
    outSAMmode_command = params.outSAMmode ? "--outSAMmode ${params.outSAMmode}" : "--outSAMmode Full"
    outSAMattributes_command = params.outSAMattributes ? "--outSAMattributes ${params.outSAMattributes}" : "--outSAMattributes All"
    outBAMcompression_command = params.outBAMcompression ? "--outBAMcompression ${params.outBAMcompression}" : "--outBAMcompression 1"
    quantMode_command = params.quantMode ? "--quantMode GeneCounts" : ""
    outFilterMismatchNoverLmax_command = params.outFilterMismatchNoverLmax ? "--outFilterMismatchNoverLmax ${params.outFilterMismatchNoverLmax}" : "--outFilterMismatchNoverLmax 0.3"
    """
    #!/bin/bash
    
    source activate rnaseq

    STAR --runMode alignReads \\
        --genomeDir ${star_idx} \\
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
        ${outFilterMismatchNmax_command} \\
        ${outFilterScoreMin_command} \\
        ${outSAMunmapped_command} \\
        ${outFilterMultimapScoreRange_command} \\
        ${alignEndsType_command} \\
        ${genomeLoad_command} \\
        ${outReadsUnmapped_command} \\
        ${outSAMmode_command} \\
        ${outSAMattributes_command} \\
        ${outBAMcompression_command} \\
        ${quantMode_command} \\
        ${outFilterMismatchNoverLmax_command} \\
        --outFileNamePrefix ${sampleID}.

    samtools index -@ 12 ${sampleID}.Aligned.sortedByCoord.out.bam
    """
}

process STAR_INDEX_RRNA {
    
    tag "${genome}_STAR_rRNA_index"

    label 'low'

    //publishDir "$params.star_index_path", mode : 'copy'
    
    input : 
        val genome
        val genome_dir
    
    output : 
        path("*")
        val("${genome_dir}"), emit : star_idx_ch

    script : 
    """
    #!/bin/bash

    source activate rnaseq
    
    STAR --runMode genomeGenerate \\
         --genomeFastaFiles ${genome} \\
         --runThreadN ${task.cpus} \\
         --genomeDir ${genome_dir}
    """

}

process STAR_REMOVE_RRNA {

    tag "${sampleID}_STAR_rRNA_remove"

    label 'medium'

    publishDir "$params.results/star_rRNA_filtered/fastq", mode : 'copy', pattern : '*.fastq.gz'
    publishDir "$params.results/star_rRNA_filtered/logs", mode : 'copy', pattern : '*Log.final.out'

    input:
        val star_idx
        tuple val(sampleID), val(fastq)   // The sample ID and the input FASTQ files (reads)

    output:
        tuple val(sampleID), path("${sampleID}_filtered.fastq.gz"), emit : filtered_reads_ch
        path("*Log.final.out")

    script:
    fastq_command = params.single_end ? "--readFilesIn ${fastq}" : "--readFilesIn ${fastq[0]} ${fastq[1]}"
    readFilesCommand = params.readFilesCommand ? "--readFilesCommand ${params.readFilesCommand}" : "zcat"

    """
    #!/bin/bash

    source activate rnaseq

    # Align reads to the rRNA reference using STAR and extract unmapped reads
    STAR --runMode alignReads \\
        --genomeDir ${star_idx} \\
        --runThreadN ${task.cpus} \\
        ${fastq_command} \\
        ${readFilesCommand} \\
        --outFileNamePrefix ${sampleID}_rRNA_removed. \\
        --outReadsUnmapped Fastx \\
        --outSAMtype None \\
        --outFilterMultimapNmax 20

    # Move unmapped reads to the filtered FASTQ file
    mv ${sampleID}_rRNA_removed.Unmapped.out.mate1 ${sampleID}_filtered.fastq.gz
    """
}

/*
    # for alignment for sailor purposes
    #STAR --runMode alignReads \\
    #    --genomeDir ${star_idx} \\
    #    --runThreadN ${task.cpus} \\
    #    ${fastq_command} \\
    #    --readFilesCommand zcat \\
    #    --alignIntronMax 10000 \\
    #    --outSAMtype BAM SortedByCoordinate \\
    #    --outFilterType BySJout \\
    #    --outFilterMultimapNmax 1 \\
    #    --outFilterIntronMotifs RemoveNoncanonical \\
    #    --outSAMprimaryFlag OneBestScore \\
    #    --outFilterScoreMinOverLread 0.0 \\
    #    --outFilterMatchNminOverLread 0.0 \\
    #    --outFilterMismatchNmax 999 \\
    #    --outSAMunmapped Within \\
    #    --outFilterMultimapScoreRange 1 \\
    #    --alignEndsType Local \\
    #    --genomeLoad NoSharedMemory \\
    #    --outReadsUnmapped Fastx \\
    #    --outSAMmode Full \\
    #    --outSAMattributes All \\
    #    --outBAMcompression 10 \\
    #    --outFileNamePrefix ${sampleID}.
*/
