#!/usr/bin/env nextflow

/*
========================================================================================
                         whole genome sequencing pipeline
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
Check bwa index is made
////////////////////////////////////////////////////////////////////
*/
genome_fasta = file("${params.genome}")
genome_name = "${genome_fasta.baseName}"
params.bwa_index_path = "${params.index}/bwa/${genome_name}"
params.bwa_index = "${params.bwa_index_path}/${genome_name}"
params.bwa_chrom_sizes = "${params.bwa_index_path}/chrom_sizes.txt"

bwa_exists = file(params.bwa_index_path).exists()

if (bwa_exists == true){
    params.bwa_build = false
} else {
    params.bwa_build = true
}

/*
////////////////////////////////////////////////////////////////////
Validate inputs
////////////////////////////////////////////////////////////////////
*/
if (params.design)    { ch_design = file(params.design, checkIfExists: true) } else { exit 1, 'Design file not specified!' }
if (params.genome)    { ch_genome = file(params.genome, checkIfExists: true) } else { exit 1, 'Genome fasta not specified!' }
if (params.results)   { ; } else { exit 1, 'Results path not specified!' }
if (params.outprefix) { ; } else { 'Outprefix not specified! Defaulting to CHIP_analysis'; params.outprefix = "CHIP_analysis" }

/*
////////////////////////////////////////////////////////////////////
Enable dls2 language --> import modules
////////////////////////////////////////////////////////////////////
*/
nextflow.enable.dsl=2

include { MRNA_SAMPLES_SHEET } from '../../modules/general/main.nf'
//include { TRIM_GALORE } from '../../modules/trimgalore/main.nf'
include { ATRIA } from '../../modules/atria/main.nf'
include { FASTQC } from '../../modules/fastqc/main.nf'
include { BWA_INDEX } from '../../modules/bwa/main.nf'
include { BWA_MEM } from '../../modules/bwa/main.nf'
include { FILTER_BAM } from '../../modules/samtools/main.nf'
include { PICARD_METRICS } from '../../modules/picard/main.nf'
include { BAM_TO_BW } from '../../modules/deeptools/main.nf'
include { DEEPVARIANT_CALL_VARIANTS } from '../../modules/deepvariant/main.nf'

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
    ATRIA( reads_ch )

    fq_ch = ATRIA.out.fq_ch

    FASTQC( reads_ch )

    /*
     * Build bwa index
     */
    if (params.bwa_build == true){
        BWA_INDEX( params.genome )
        bwa_index_ch = BWA_INDEX.out.bwa_index_ch
    } else {
        bwa_index_ch = params.bwa_index 
    }

    /*
     * Align with bwa
     */
    BWA_MEM( bwa_index_ch, fq_ch )

    /*
     * Filter merged bam files
     * Some code used and inspired from https://github.com/nf-core/chipseq
     */

    if ( params.filter_bam ){
        bam_ch = BWA_MEM.out.bwa_bam_ch
        FILTER_BAM( bam_ch )
        processed_bam_ch = FILTER_BAM.out.filtered_bams_ch
        bam_to_bw_input_ch = FILTER_BAM.out.bam_to_bw_input
    } else {
        processed_bam_ch = BWA_MEM.out.bwa_bam_ch
        bam_to_bw_input_ch = BWA_MEM.out.bwa_bam_bai_ch
    }

    /*
     * BAM --> BW
     */
    BAM_TO_BW( bam_to_bw_input_ch )


    /* 
     * call variants 
     */
    // DeepVariant

    DEEPVARIANT_CALL_VARIANTS( bam_to_bw_input_ch, params.genome )

    // Repeat Expansion
    // Indel caller
    
}
