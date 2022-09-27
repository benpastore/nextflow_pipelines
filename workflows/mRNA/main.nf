#!/usr/bin/env nextflow

/*
 * Set bool for single/paired end library
 */
if (params.library == 'single'){
    params.single_end = true
} else {
    params.single_end = false
}   

/*
 * Set command for gene_ids
 */
if (params.annotation.gene_ids == false) {
    params.annotation.gene_ids = ""
} else {
    params.gene_ids = "-a ${params.annotation.gene_ids}"
}

/*
 * Set path to bin
 */
params.bin = "${baseDir}/bin"

/*
 * Check if star genome exists
 */
genome_fasta = file("${params.annotation.genome}")
genome_name = "${genome_fasta.baseName}"

params.star_index = "${baseDir}/index/star/${genome_name}"
star_path = file("${params.star_index}")
star_exists = "${star_path.exists()}"

if (star_exists == true){
    params.star_build = false
} else {
    params.star_build = true
}

/*
 * Check if salmon index exists
 */
if (params.salmon.run == true){
    transcript_fasta = file("${params.salmon.reference}")
    transcript_name = "${transcript_fasta.baseName}"
    
    params.salmon_index = "${baseDir}/index/salmon/${transcript_name}"
    salmon_path = file("${params.salmon_index}")
    salmon_exists = "${salmon_path.exists()}"
    
    if (salmon_exists == true){
        params.salmon_build = false
    } else {
        params.salmon_build = true
    }
} else {
    params.salmon_build = false
}

/*
 * Enable dsl2 syntax
 */
nextflow.enable.dsl=2

/*
 * Import modules
 */
include {
    TRIM_GALORE;
    FASTQC;
    STAR_INDEX;
    SALMON_INDEX;
    STAR_ALIGN;
    SALMON;
    BAM_TO_BW;
    COUNTER;
    RBIND_COUNTS;
    MASTER_TABLE;
    MULTI_QC
} from './modules.nf' 

/*
 * Pipeline logic
 */
workflow {

    /*
     * Set reads channel
     */
    if (params.single_end == true){
        reads = Channel.fromPath( params.reads ).ifEmpty{ error "No reads in reads directory" }.map { file -> tuple(file.simpleName, file) }
    } else {
        reads = Channel.fromFilePairs( params.reads ).ifEmpty{ error "No reads in reads directory" }
    }

    /*
     * Trim reads of adapters and low quality sequences
     */
    TRIM_GALORE( reads )
    TRIM_GALORE.out.fq_ch.tap{ fq_ch }
    
    /*
     * Preform FASTQC on raw reads
     */
    FASTQC( reads )

    /*
     * STAR Alignment.
     * 1. Build star index if necessary
     * 2. Preform alignment
     */
    if (params.star_build == true){
        STAR_INDEX ( params.annotation.genome, params.annotation.gtf)
        STAR_INDEX.out.star_idx_ch.tap{ star_idx_ch }
    } else {
        star_idx_ch = star_path
    }
    
    STAR_ALIGN(star_idx_ch, params.annotation.gtf, fq_ch)
    STAR_ALIGN.out.star_bams_ch.tap{ bams }

    /*
     * SALMON Alignment
     * 1. Only run if indicated to run
     * 2. Build samlmon index if necessary
     * 3. Preform salmon mapping 
     */
    if (params.salmon.run == true){
        if (params.salmon_build == true){
            SALMON_INDEX ( params.salmon.reference )
            SALMON_INDEX.out.salmon_idx_ch.tap{ salmon_idx_ch }
        } else {
            salmon_idx_ch = salmon_path
        }
    }

    if (params.salmon.run == true){
        SALMON( salmon_idx_ch, fq_ch )
        SALMON.out.salmon_quant_ch.tap{ salmon_quant_ch }
    }

    /*
     * BAM --> BW using deeptools bam coverage
     */
    BAM_TO_BW( bams )

    /*
     * Feature Counts
     */
    COUNTER( bams, params.annotation.gtf )
    COUNTER.out.feature_counts_ch.tap{ feature_counts_ch }
    
    /*
     * Join counts table from salmon + feature counts, if neccesary 
     */
    if (params.salmon.run == true) {
        feature_counts_ch.join(salmon_quant_ch).tap{ counts_ch }
        RBIND_COUNTS( counts_ch )
        RBIND_COUNTS.out.counts_ch.tap{ master_table_input }
    } else {
        master_table_input = feature_counts_ch
    }

    /*
     * Make master table
     */
    MASTER_TABLE(
        params.project.name,
        master_table_input.collect(),
        params.gene_ids
    )
    
    /*
     * MultiQC
     */
    MULTI_QC(MASTER_TABLE.out.master_tables_ch, params.results )

}