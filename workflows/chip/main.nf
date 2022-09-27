#!/usr/bin/env nextflow

/*
========================================================================================
                         CHIP sequencing pipeline
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

bwa_exists = file(params.bwa_index_path).exists()

if (bwa_exists == true){
    params.bwa_build = false
} else {
    params.bwa_build = true
}

/*
////////////////////////////////////////////////////////////////////
Enable dls2 language --> import modules
////////////////////////////////////////////////////////////////////
*/
nextflow.enable.dsl=2

include { CHIP_SAMPLES_SHEET } from '../../modules/chip_samples_sheet/main.nf'
include { TRIM_GALORE } from '../../modules/trimgalore/main.nf'
include { FASTQC } from '../../modules/fastqc/main.nf'
include { BWA_INDEX } from '../../modules/bwa/main.nf'
include { BWA_MEM } from '../../modules/bwa/main.nf'
include { MERGE_BAM } from '../../modules/picard/main.nf'
include { FILTER_MERGE_BAM } from '../../modules/picard/main.nf'
include { PICARD_METRICS } from '../../modules/picard/main.nf'
include { MACS2 } from '../../modules/macs2/main.nf'
include { BAM_TO_BW } from '../../modules/deeptools/main.nf'
include { BAM_COMPARE } from '../../modules/deeptools/main.nf'
include { MULTI_QC } from '../../modules/multiqc/main.nf'

/*
////////////////////////////////////////////////////////////////////
Workflow
////////////////////////////////////////////////////////////////////
*/
workflow {

    /*
     * Process CHIP Samples sheet input
     */
    CHIP_SAMPLES_SHEET( params.design )

    if (params.single_end == false ){
        CHIP_SAMPLES_SHEET.out.fq_ch
            .splitCsv(header:true, sep:',')
            .map { row -> [ row.simple_name, [ file(row.R1,checkIfExists: true), file(row.R2,checkIfExists: true) ] ] }
            .tap{ reads_ch }
    } else {
        CHIP_SAMPLES_SHEET.out.fq_ch
            .splitCsv(header:true, sep:',')
            .map { row -> [ row.simple_name, file(row.R1,checkIfExists: true) ] }
            .tap{ reads_ch }
    }
    
    CHIP_SAMPLES_SHEET.out.control_ch
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.group, row.control_group ] }
        .tap{ control_ch }

    CHIP_SAMPLES_SHEET.out.replicates_ch
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.simple_name, row.group ] }
        .tap{ replicates_ch }

    /*
     * Trim reads of adapters and low quality sequences
     */
    TRIM_GALORE( reads_ch )
    TRIM_GALORE.out.fq_ch.tap{ fq_ch }

    FASTQC( reads_ch )

    /*
     * Build bwa index
     */
    if (params.bwa_build == true){
        BWA_INDEX( params.genome )
        BWA_INDEX.out.tap{ bwa_index_ch }
    } else {
        bwa_index_ch = params.bwa_index 
    }

    /*
     * Align with bwa
     */
    BWA_MEM( bwa_index_ch, fq_ch )
    BWA_MEM.out.bwa_bam_ch.tap{ bam_ch }

    bam_ch
        .join(replicates_ch)
        .groupTuple(by:2)
        .tap{ merge_bams }

    /*
     * Merge BAM files 
     * Some code used and inspired from https://github.com/nf-core/chipseq
     */
    MERGE_BAM( merge_bams )

    /*
     * Filter merged bam files
     * Some code used and inspired from https://github.com/nf-core/chipseq
     */
    FILTER_MERGE_BAM( MERGE_BAM.out.merged_bams_ch )
    FILTER_MERGE_BAM.out.filtered_bams_ch.tap{ filtered_bams_ch }

    filtered_bams_ch
        .join(control_ch)
        .map { it -> it[2,0,1] }
        .filter { it[0] != it[1] }
        .combine(filtered_bams_ch)
        .filter { it[0] == it[3] }
        .map { it -> it[1,2,0,4]}
        .tap { bam_compare_ch }

    /*
     * BAM --> BW
     */
    BAM_TO_BW( FILTER_MERGE_BAM.out.bam_to_bw_input )

    /*
     * MACS2
     */
    MACS2( bam_compare_ch ) 

    /*
     * BAM Compare
     */
    BAM_COMPARE( bam_compare_ch )

    /*
     * Picard metrics
     * Some code used and inspired from https://github.com/nf-core/chipseq
     */
    PICARD_METRICS(filtered_bams_ch, params.genome)
    
}
