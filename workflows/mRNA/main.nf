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

/*
////////////////////////////////////////////////////////////////////
Check star index is made
////////////////////////////////////////////////////////////////////
*/
genome_fasta = file("${params.genome}")
genome_name = "${genome_fasta.baseName}"
params.star_index_path = "${params.index}/star/${genome_name}"
params.star_index = "${params.star_index_path}"
params.star_chrom_sizes = "${params.star_index_path}/chrNameLength.txt"

star_exists = file(params.star_index_path).exists()

if (star_exists == true){
    params.star_build = false
} else {
    params.star_build = true
}

/*
////////////////////////////////////////////////////////////////////
Check salmon index is made
////////////////////////////////////////////////////////////////////
*/
if (params.salmon) {

    genome_fasta = file("${params.salmon}")
    genome_name = "${genome_fasta.baseName}"
    params.salmon_index_path = "${params.index}/salmon/${genome_name}"
    params.salmon_index = "${params.index}/salmon"

    salmon_exists = file(params.salmon_index_path).exists()

    if (salmon_exists == true){
        params.salmon_build = false
    } else {
        params.salmon_build = true
    }

}

/*
////////////////////////////////////////////////////////////////////
Validate inputs
////////////////////////////////////////////////////////////////////
*/
if (params.design)          { ch_design = file(params.design, checkIfExists: true) } else { exit 1, 'Design file not specified!' }
if (params.genome)          { ch_genome = file(params.genome, checkIfExists: true) } else { exit 1, 'Genome fasta not specified!' }
if (params.results)         { ; } else { exit 1, 'Results path not specified!' }
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
Subworkflow
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


/*
////////////////////////////////////////////////////////////////////
Workflow
////////////////////////////////////////////////////////////////////
*/
workflow {

    /*
     * Trim reads of adapters and low quality sequences
     */
    MRNA_SAMPLES_SHEET( params.design )

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

    /*
     * Trim reads of adapters and low quality sequences
     */
    TRIM_GALORE( reads_ch )
    fq_ch = TRIM_GALORE.out.fq_ch
    
    /*
     * Preform FASTQC on raw reads
     */
    FASTQC( reads_ch )

    /*
     * STAR Alignment.
     * 1. Build star index if necessary
     * 2. Preform alignment
     */
    if ( params.star_build ){
        STAR_INDEX ( params.genome, params.gtf, params.star_index_path)
        star_idx_ch = STAR_INDEX.out.star_idx_ch
    } else {
        star_idx_ch = params.star_index_path
    }
    
    STAR_ALIGN(star_idx_ch, params.gtf, fq_ch)

    bam_ch = STAR_ALIGN.out.star_bams_ch

    /*
     * Filter STAR alignment
     */
    
    if ( params.filter_bam ){
        FILTER_BAM( bam_ch )
        processed_bam_ch = FILTER_BAM.out.filtered_bams_ch
        bam_to_bw_input_ch = FILTER_BAM.out.bam_to_bw_input
    } else {
        processed_bam_ch = bam_ch
        bam_to_bw_input_ch = bam_ch
    }

    /*
     * BAM --> BW using deeptools bam coverage
     */
    BAM_TO_BW( bam_to_bw_input_ch )

    /*
     * Merge BW files for a given condition
     */
    if ( params.merge_bw ){

        merge_bws = BAM_TO_BW
            .out
            .bws_ch
            .join(replicates_ch)
            .groupTuple(by:2)

        MERGE_BW(merge_bws, params.star_chrom_sizes)
    
    }

    /*
     * SALMON Alignment
     * 1. Only run if indicated to run
     * 2. Build samlmon index if necessary
     * 3. Preform salmon mapping 
     */
    if ( params.salmon ){
        if ( params.salmon_build ){
            SALMON_INDEX ( params.salmon )
            salmon_idx_ch = SALMON_INDEX.out.salmon_idx_ch
        } else {
            salmon_idx_ch = params.salmon_index_path
        }

        SALMON( salmon_idx_ch, fq_ch )
        salmon_quant_ch = SALMON.out.salmon_quant_ch
    }

    /*
     * Feature Counts
     */
    if ( params.feature_counts ) {
        FEATURECOUNTS( processed_bam_ch, params.gtf )
        feature_counts_ch = FEATURECOUNTS.out.feature_counts_ch
        
        /*
        * Join counts table from salmon + feature counts, if neccesary 
        */
        if (params.salmon) {
            counts_ch = feature_counts_ch.join(salmon_quant_ch)

            RBIND_COUNTS( counts_ch )
            
            master_table_input = RBIND_COUNTS.out.counts_ch
        
        } else {
            master_table_input = feature_counts_ch
        
        }

        params.project_name = params.outprefix ? "${params.outprefix}" : "mRNA_analysis"
        /*
        * Make master table
        */
        MASTER_TABLE(
            params.project_name,
            master_table_input.collect()
        )
    }

}


