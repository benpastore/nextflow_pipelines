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

include { CHIP_SAMPLES_SHEET } from '../../modules/general/main.nf'
include { TRIM_GALORE } from '../../modules/trimgalore/main.nf'
include { FASTQC } from '../../modules/fastqc/main.nf'
include { BWA_INDEX } from '../../modules/bwa/main.nf'
include { BWA_MEM } from '../../modules/bwa/main.nf'
include { MERGE_BAM } from '../../modules/picard/main.nf'
include { FILTER_BAM } from '../../modules/picard/main.nf'
include { PICARD_METRICS } from '../../modules/picard/main.nf'
include { MACS2 } from '../../modules/macs2/main.nf'
include { BAM_TO_BW } from '../../modules/deeptools/main.nf'
include { BAM_COMPARE } from '../../modules/deeptools/main.nf'
include { MULTI_QC } from '../../modules/multiqc/main.nf'
include { CHIP_INTERSECT } from '../../modules/bedtools/main.nf'

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
        reads_ch = CHIP_SAMPLES_SHEET.out.fq_ch
            .splitCsv(header:true, sep:',')
            .map { row -> [ row.simple_name, [ file(row.R1,checkIfExists: true), file(row.R2,checkIfExists: true) ] ] }
            
    } else {
        reads_ch = CHIP_SAMPLES_SHEET.out.fq_ch
            .splitCsv(header:true, sep:',')
            .map { row -> [ row.simple_name, file(row.R1,checkIfExists: true) ] }
            
    }
    
    control_ch = CHIP_SAMPLES_SHEET.out.control_ch
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.group, row.control_group ] }
        

    replicates_ch = CHIP_SAMPLES_SHEET.out.replicates_ch
        .splitCsv(header:true, sep:',')
        .map { row -> [ row.simple_name, row.group ] }
        

    /*
     * Trim reads of adapters and low quality sequences
     */
    TRIM_GALORE( reads_ch )
    fq_ch = TRIM_GALORE.out.fq_ch

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
    bam_ch = BWA_MEM.out.bwa_bam_ch

    merge_bams = bam_ch
        .join(replicates_ch)
        .groupTuple(by:2)

    /*
     * Merge BAM files 
     * Some code used and inspired from https://github.com/nf-core/chipseq
     */
    MERGE_BAM( merge_bams )

    /*
     * Filter merged bam files
     * Some code used and inspired from https://github.com/nf-core/chipseq
     */
    FILTER_BAM( MERGE_BAM.out.merged_bams_ch )
    filtered_bams_ch = FILTER_BAM.out.filtered_bams_ch

    bam_compare_ch = filtered_bams_ch
        .join(control_ch)
        .map { it -> it[2,0,1] }
        .filter { it[0] != it[1] }
        .combine(filtered_bams_ch)
        .filter { it[0] == it[3] }
        .map { it -> it[1,2,0,4]}

    /*
     * BAM --> BW
     */
    BAM_TO_BW( FILTER_BAM.out.bam_to_bw_input )

    /*
     * MACS2
     */
    MACS2( bam_compare_ch ) 

    /*
     * BAM Compare
     */
    BAM_COMPARE( bam_compare_ch )

    /*
     * intersect bam to peaks
     */
    peaks_ch = MACS2.out.macs2_peaks_ch
    
    peak_intersect_ch = filtered_bams_ch
        .join(peaks_ch)
    
    CHIP_INTERSECT( peak_intersect_ch )

    /*
     * Picard metrics
     * Some code used and inspired from https://github.com/nf-core/chipseq
     */
    PICARD_METRICS(filtered_bams_ch, params.genome)
    
}
