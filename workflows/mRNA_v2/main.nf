#!/usr/bin/env nextflow

/*
========================================================================================
                         mRNA sequencing pipeline
========================================================================================
Ben Pastore
pastore.28@osu.edu
----------------------------------------------------------------------------------------
*/

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:

    nextflow main.nf -profile cluster --design design.csv 

    Mandatory arguments:
    --design [file]                     Comma-separated file containing information about the samples in the experiment (see docs/usage.md) (Default: './design.csv')
    --results [file]                    Path to results directory
    --single_end [bool]                 IF using single end sequencing must specify as true.
    --outprefix [name]                  outprefix to be assigned to output files when necessary, default is "mRNA_analysis"
    -profile [str]                      Configuration profile to use. Can use local / cluster
    """.stripIndent()
}

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

params.bin = "${params.base}/../../bin"
params.index = "${params.base}/../../index"
params.date = new Date().format( 'yyyyMMdd' )
params.results = "${params.outdir}/${params.date}"

/*
////////////////////////////////////////////////////////////////////
Validate inputs
////////////////////////////////////////////////////////////////////
*/
if (params.design)          { ch_design = file(params.design, checkIfExists: true) } else { exit 1, 'Design file not specified!' }
if (params.genome)          { ch_genome = file(params.genome, checkIfExists: true) } else { exit 1, 'Genome fasta not specified!' }
if (params.outprefix)       { ; } else { 'Outprefix not specified! Defaulting to mRNA_analysis' }

/*
////////////////////////////////////////////////////////////////////
set project name
////////////////////////////////////////////////////////////////////
*/
params.project_name = params.outprefix ? "${params.outprefix}" : "mRNA_analysis"

/*
////////////////////////////////////////////////////////////////////
Enable dls2 language --> import modules
////////////////////////////////////////////////////////////////////
*/
nextflow.enable.dsl=2

include { MRNA_SAMPLES_SHEET } from '../../modules/general/main.nf'
include { TRIM_GALORE } from '../../modules/trimgalore/main.nf'
include { FASTQC } from '../../modules/fastqc/main.nf'
include { STAR_INDEX } from '../../modules/star/main.nf'
include { SALMON_INDEX } from '../../modules/salmon/main.nf'
include { SALMON } from '../../modules/salmon/main.nf'
include { STAR_ALIGN } from '../../modules/star/main.nf'
include { PICARD_METRICS } from '../../modules/picard/main.nf'
include { BAM_TO_BW } from '../../modules/deeptools/main.nf'
include { MERGE_BW } from '../../modules/deeptools/main.nf'
include { MULTI_QC } from '../../modules/multiqc/main.nf'
include { MERGE_BAM } from '../../modules/picard/main.nf'
include { FILTER_BAM } from '../../modules/samtools/main.nf'
include { FEATURECOUNTS } from '../../modules/featurecounts/main.nf'
include { RBIND_COUNTS } from '../../modules/general/main.nf'
include { MASTER_TABLE } from '../../modules/general/main.nf'

/*
////////////////////////////////////////////////////////////////////
Import subworkflows ???
////////////////////////////////////////////////////////////////////
*/

/*
////////////////////////////////////////////////////////////////////
Subworkflows
////////////////////////////////////////////////////////////////////
*/
workflow FEATURE_COUNT {
    /*
     * Parse bam files into channel taking bam file path and simple name
     */
    bams_ch = Channel
        .fromPath( "$params.bam_path/*.bam" )
        .map { fastq -> [ fastq.SimpleName, fastq] }
    
    bams_ch.view()

    /*
     * Feature Counts
     */
    FEATURECOUNTS( bams_ch, params.gtf )
    
    feature_counts_ch = FEATURECOUNTS.out.feature_counts_ch
    
    master_table_input = feature_counts_ch
    
    /*
     * set project name
     */
    params.project_name = params.outprefix ? "${params.outprefix}" : "feature_counting"

    /*
     * Make master table
     */
    MASTER_TABLE(
        params.project_name,
        master_table_input.collect()
    )
}

workflow read_data_fq {

    take :
        data

    main : 
        MRNA_SAMPLES_SHEET( data )

        if (params.single_end == false ){
            reads_ch = MRNA_SAMPLES_SHEET.out.fq_ch
                .splitCsv(header:true, sep:',')
                .map { row -> [ row.simple_name, [ file(row.R1,checkIfExists: true), file(row.R2,checkIfExists: true) ] ] }
                
        } else {
            reads_ch = MRNA_SAMPLES_SHEET.out.fq_ch
                .splitCsv(header:true, sep:',')
                .map { row -> [ row.simple_name, file(row.R1,checkIfExists: true) ] }
                
        }
        
        replicates_ch = MRNA_SAMPLES_SHEET.out.replicates_ch
            .splitCsv(header:true, sep:',')
            .map { row -> [ row.simple_name, row.group ] }
    
    emit : 
        reads = reads_ch
        replicates = replicates_ch
}

workflow trim {

    take : 
        data

    main :
        /*
        * Trim reads of adapters and low quality sequences
        */
        TRIM_GALORE( data )

        fq_ch = TRIM_GALORE.out.fq_ch

        if ( params.fastqc ){
            FASTQC( data )
        }
    
    emit : 
        fqs = fq_ch
}

workflow star_align {

    take : 
        data 
    
    main : 
        genome_fasta = file("${params.genome}")
        genome_name = "${genome_fasta.baseName}"
        star_index_path = "${params.index}/star/${genome_name}"
        star_index = "${star_index_path}"
        star_chrom_sizes = "${star_index_path}/chrNameLength.txt"

        star_exists = file(star_chrom_sizes).exists()
        println star_exists
        if (star_exists == true){
            star_build = false
        } else {
            star_build = true
        }

        if (star_build == true){

            STAR_INDEX( params.genome, params.gtf, star_index_path )
            star_index_ch = STAR_INDEX.out.star_idx_ch

        } else {

            star_index_ch = star_index 
        
        }

        STAR_ALIGN( star_index_ch, params.gtf, data )
    
    emit : 
        bams = STAR_ALIGN.out.star_bam_bai_ch
        chrom_sizes = "${star_index_path}/chrNameLength.txt"

}

workflow filter_bam {

    take : data

    main : 

        FILTER_BAM( data )

        processed_bam_ch = FILTER_BAM.out.filtered_bams_ch
        
        bam_to_bw_input_ch = FILTER_BAM.out.bam_to_bw_input

    emit : 
        filter_bam_bai = bam_to_bw_input_ch

}

workflow bam_to_bw {

    take : data

    main : 
        BAM_TO_BW( data )
    
    emit : 
        bws = BAM_TO_BW.out.bws_ch

}

workflow merge_bw { 

    take : 
        data
        replicates
        chrom_sizes
    
    main : 
        merge_bws = data
            .join(replicates)
            .groupTuple(by:2)

        MERGE_BW(merge_bws, chrom_sizes)

}

workflow salmon {

    take : 
        data
    
    main :
        genome_fasta = file("${params.salmon}")
        genome_name = "${genome_fasta.baseName}"
        salmon_index_path = "${params.index}/salmon/${genome_name}"
        salmon_index = "${params.index}/salmon"

        salmon_exists = file(salmon_index_path).exists()
        if (salmon_exists == true){
            salmon_build = false
        } else {
            salmon_build = true
        }

        if ( salmon_build ){
            SALMON_INDEX ( params.salmon, salmon_index_path )
            salmon_idx_ch = SALMON_INDEX.out.salmon_idx_ch
        } else {
            salmon_idx_ch = salmon_index_path
        }

        SALMON( salmon_idx_ch, data )
        salmon_quant_ch = SALMON.out.salmon_quant_ch
    
    emit : 
        salmon_quant = SALMON.out.salmon_quant_ch

}

workflow feature_count {

    take : 
        data 
    
    main : 
        FEATURECOUNTS( data, params.gtf )
    
    emit : 
        counts = FEATURECOUNTS.out.feature_counts_ch

}

workflow rbind_counts { 
    take : 
        data
    main : 
        RBIND_COUNTS( data )
    emit : 
        counts = RBIND_COUNTS.out.counts_ch
}

workflow master_table { 

    take : 
        name
        data
    
    main : 
        MASTER_TABLE(
                params.project_name,
                data.collect()
            )

}

/*
////////////////////////////////////////////////////////////////////
Workflow
////////////////////////////////////////////////////////////////////
*/
workflow {

    // read fastq data
    read_data_fq( params.design )

    // trim 
    trim( read_data_fq.out.reads )

    // star align
    star_align( trim.out.fqs )
    
    // filter bam
    if ( params.filter_bam ){
        filter_bam( star_align.out.bams )
        bams = filter_bam.out.filter_bam_bai
    } else {
        bams = star_align.out.bams
    }

    // bam to bigwig
    bam_to_bw( bams )

    // merge bams
    if ( params.merge_bw ) {
        merge_bw( bam_to_bw.out.bws, read_data_fq.out.replicates, star_align.out.chrom_sizes )
    }

    /*
     * SALMON mapping
     */
    if ( params.salmon ) {
        salmon( trim.out.fqs,  )
    }

    /*
     * Feature Counts
     */
    if ( params.feature_counts ) {
        
        feature_count( bams )
        /*
        * Join counts table from salmon + feature counts, if neccesary 
        */
        if (params.salmon) {
            counts_ch = feature_count.out.counts.join(salmon.out.salmon_quant)
            rbind_counts(counts_ch)
            master_table_input = rbind_counts.out.counts
        } else {
            master_table_input = feature_count.out.counts
        }

        params.project_name = params.outprefix ? "${params.outprefix}" : "mRNA_analysis"
        /*
        * Make master table
        */
        master_table( params.project_name, master_table_input )
    }

}


