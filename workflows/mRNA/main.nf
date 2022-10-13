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
    params.salmon_index = "${params.salmon_index_path}"

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
if (params.design)    { ch_design = file(params.design, checkIfExists: true) } else { exit 1, 'Design file not specified!' }
if (params.genome)    { ch_genome = file(params.genome, checkIfExists: true) } else { exit 1, 'Genome fasta not specified!' }
if (params.results)   { ; } else { exit 1, 'Results path not specified!' }
if (params.outprefix) { ; } else { 'Outprefix not specified! Defaulting to mRNA_analysis' }

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
include { MULTI_QC } from '../../modules/multiqc/main.nf'
include { MERGE_BAM } from '../../modules/picard/main.nf'
include { FILTER_BAM } from '../../modules/picard/main.nf'
include { FEATURECOUNTS } from '../../modules/featurecounts/main.nf'
include { RBIND_COUNTS } from '../../modules/general/main.nf'
include { MASTER_TABLE } from '../../modules/general/main.nf'

/*
////////////////////////////////////////////////////////////////////
Workflow
////////////////////////////////////////////////////////////////////
*/
workflow {

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
        STAR_INDEX ( params.genome, params.gtf)
        star_idx_ch = STAR_INDEX.out.star_idx_ch
    } else {
        star_idx_ch = params.star_index_path
    }
    
    STAR_ALIGN(star_idx_ch, params.gtf, fq_ch)

    bam_ch = STAR_ALIGN.out.star_bams_ch

    /*
     * Filter STAR alignment
     */
    FILTER_BAM( bam_ch )

    filter_bam_ch = FILTER_BAM.out.filtered_bams_ch

    merge_bams = filter_bam_ch
        .join(replicates_ch)
        .groupTuple(by:2)

    /*
     * Merge BAM files 
     * Some code used and inspired from https://github.com/nf-core/chipseq
     */
    MERGE_BAM( merge_bams )

    /*
     * BAM --> BW using deeptools bam coverage
     */
    merged_bams = MERGE_BAM.out.merged_bams_ch
    BAM_TO_BW( merged_bams )

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
    FEATURECOUNTS( filter_bam_ch, params.gtf )
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
        master_table_input.collect(),
    )

}


