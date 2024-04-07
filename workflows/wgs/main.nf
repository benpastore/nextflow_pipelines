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
params.snpEff_data = "${params.bin}/snpEff"
/*
////////////////////////////////////////////////////////////////////
Check bwa index is made
////////////////////////////////////////////////////////////////////
*/
/*
* This part is setting up for the default values of the input (we can always change these value in the command line 
* in the run.sh in the sam_cendr folder). However, we DO NOT want to change these params
*/

/*
////////////////////////////////////////////////////////////////////
Validate inputs (these are things that being specified in run.sh in projects/2024_CeNDR)
////////////////////////////////////////////////////////////////////
*/
//if (params.design)    { ch_design = file(params.design, checkIfExists: true) } else { exit 1, 'Design file not specified!' }
//if (params.genome)    { ch_genome = file(params.genome, checkIfExists: true) } else { exit 1, 'Genome fasta not specified!' }
//if (params.results)   { ; } else { exit 1, 'Results path not specified!' }
//if (params.outprefix) { ; } else { 'Outprefix not specified! Defaulting to CHIP_analysis'; params.outprefix = "CHIP_analysis" }

genome_fasta = file("${params.genome}")
params.genome_name = "${genome_fasta.baseName}"
params.index_path = "${params.index}/${params.aligner}/${params.genome_name}"

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
include { INDEX_BAM } from '../../modules/samtools/main.nf'
include { PICARD_METRICS } from '../../modules/picard/main.nf'
include { BAM_TO_BW } from '../../modules/deeptools/main.nf'
include { DEEPVARIANT_CALL_VARIANTS } from '../../modules/deepvariant/main.nf'
include { EXPANSION_HUNTER } from '../../modules/expansion_hunter/main.nf'
include { GET_CONTIGS } from '../../modules/gatk/main.nf'
include { GATK_PREPARE_GENOME } from '../../modules/gatk/main.nf'
include { GATK_PROCESS_BAM } from '../../modules/gatk/main.nf'
include { GATK_CALL_VARIANTS } from '../../modules/gatk/main.nf'
include { CONCAT_STRAIN_GVCFS } from '../../modules/gatk/main.nf'
include { MAKE_SAMPLE_MAP } from '../../modules/gatk/main.nf'
include { IMPORT_GENOME_DB } from '../../modules/gatk/main.nf'
include { CONCATENATE_VCF } from '../../modules/gatk/main.nf'
include { GATK_GENOTYPE_COHORT } from '../../modules/gatk/main.nf'
include { MULTIQC } from '../../modules/gatk/main.nf'
include { BUILD_SNPEFF_DB } from '../../modules/snpeff/main.nf'
include { RUN_SNPEFF } from '../../modules/snpeff/main.nf'
include { SOFT_FILTER } from '../../modules/gatk/main.nf'
include { HARD_FILTER } from '../../modules/gatk/main.nf'
include { HTML_REPORT } from '../../modules/gatk/main.nf'
include { DOWNLOAD_BAM } from '../../modules/download/main.nf'
include { PREPROCESS_SNPEFF } from '../../modules/snpeff/main.nf'

/*
////////////////////////////////////////////////////////////////////
set date variable
////////////////////////////////////////////////////////////////////
*/
params.date = new Date().format( 'yyyyMMdd' )
println params.date

/*
////////////////////////////////////////////////////////////////////
Subworkflows
////////////////////////////////////////////////////////////////////
*/

workflow download_bam {

    take : 
        data
    
    main : 
        DOWNLOAD_BAM( data )

    emit :
        bams = DOWNLOAD_BAM.out.bam_ch 

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
        ATRIA( data )

        fq_ch = ATRIA.out.fq_ch

        if ( params.fastqc ){
            FASTQC( data )
        }
    
    emit : 
        fqs = fq_ch
}

workflow index_bam {

    take : data 

    main : 
        INDEX_BAM( data )
    
    emit : 
        bams = INDEX_BAM.out

}

workflow bwa_align {

    take : data

    main :
        genome_fasta = file("${params.genome}")
        genome_name = "${genome_fasta.baseName}"
        bwa_index_path = "${params.index}/bwa/${genome_name}"
        bwa_index = "${bwa_index_path}/${genome_name}"
        bwa_chrom_sizes = "${bwa_index_path}/chrom_sizes.txt"

        bwa_exists = file(bwa_index_path).exists()

        if (bwa_exists == true){
            bwa_build = false
        } else {
            bwa_build = true
        }

        if (bwa_build == true){
            BWA_INDEX( params.genome, bwa_index )
            bwa_index_ch = BWA_INDEX.out.bwa_index_ch
        } else {
            bwa_index_ch = bwa_index 
        }

        BWA_MEM( bwa_index_ch, data )
    
    emit :
        bams = BWA_MEM.out.bwa_bam_bai_ch
}

workflow star_align {

    take : 
        data 
    
    main : 
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

        if (params.star_build == true){

            STAR_INDEX( params.genome, params.gtf, params.star_index_path )
            star_index_ch = STAR_INDEX.out.star_idx_ch

        } else {

            star_index_ch = params.star_index 
        
        }

        STAR_ALIGN( star_index_ch, params.gtf, data )
    
    emit : 
        bams = STAR_ALIGN.out.star_bam_bai_ch

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

}

workflow deepvariant {

    take : data 

    main :
        DEEPVARIANT_CALL_VARIANTS( data, params.genome )

    emit : 
        vcfs = DEEPVARIANT_CALL_VARIANTS.out.vcf

}

workflow gatk {

    take : data 

    main : 
        /*
        * Get contigs
        */
        GET_CONTIGS( data.first() )

        // get contigs: I, II, III, IV, V, X
        contigs = GET_CONTIGS.out.splitText { it.strip() }

        GATK_PREPARE_GENOME( params.genome )

        GATK_PROCESS_BAM( data )

        variant_caller_input = GATK_PROCESS_BAM
            .out
            .processed_bam_ch
            .combine( contigs )
            .combine( GATK_PREPARE_GENOME.out.processed_genome )

        // GATK
        GATK_CALL_VARIANTS( variant_caller_input )

        // concat vcfs, combine chromsome calls for each strain
        CONCAT_STRAIN_GVCFS( GATK_CALL_VARIANTS.out.groupTuple() )
    
    emit : 
        vcfs = CONCAT_STRAIN_GVCFS.out.concat_vcf_ch.collect().toList()
        refs = GATK_PREPARE_GENOME.out.processed_genome
        contigs = GET_CONTIGS.out.contigs.splitText { it.strip() }
        contig_file = GET_CONTIGS.out.contigs

}

workflow get_contigs {

    take : 
        data 

    main : 
        GET_CONTIGS( data )
    
    emit : 
        contigs = GET_CONTIGS.out.contigs.splitText { it.strip() }
        contig_file = GET_CONTIGS.out.contigs

}

workflow gatk_prepare_genome {
    
    take : 
        data
    
    main : 
        GATK_PREPARE_GENOME( data )
    
    emit : 
        refs = GATK_PREPARE_GENOME.out.processed_genome
    
}

workflow gatk_genotype_cohort {

    take : 
        vcfs
        refs
        contigs 
    
    main : 
        // write sample map 

        sample_vcfs = vcfs
            .flatten()
            .map { it -> "$it" }
            .collectFile( name : "sample_vcfs.tsv", storeDir : "$params.results", newLine : true)

        MAKE_SAMPLE_MAP( sample_vcfs )

        genome_db_input = MAKE_SAMPLE_MAP
            .out
            .sample_map
            .combine( contigs )

        // make gatk db
        IMPORT_GENOME_DB( genome_db_input )

        // genotype cohort gvcf db
        IMPORT_GENOME_DB
            .out
            .genome_db_output
            .combine( refs ) | \
            GATK_GENOTYPE_COHORT

        // combine vcf
        GATK_GENOTYPE_COHORT.out.genotype_cohort_gvcf_db_out.collect() | CONCATENATE_VCF
    
    emit : 
        vcf_tbi = CONCATENATE_VCF.out.vcf_tbi
        vcfs = CONCATENATE_VCF.out.vcf

}

workflow gatk_filters {

    take : 
        vcfs 
        refs 
        contigs
    
    main : 
        vcfs.combine( refs ) | SOFT_FILTER

        SOFT_FILTER.out.soft_filter_vcf.combine(contigs) | HARD_FILTER
    
    emit : 
        vcf_filter_stats = SOFT_FILTER.out.soft_vcf_stats.concat(HARD_FILTER.out.hard_vcf_stats).collect()
        vcfs = SOFT_FILTER.out.vcf.mix(HARD_FILTER.out.hard_filter_vcf).flatten()

}

workflow snpeff {

    take : 
        vcfs 
    
    main : 

        PREPROCESS_SNPEFF( params.genome, params.gtf, params.protein, params.cds, params.organism )

        BUILD_SNPEFF_DB( PREPROCESS_SNPEFF.out.done )

        RUN_SNPEFF( vcfs, BUILD_SNPEFF_DB.out.db_name )

}

workflow expansion_hunter {

    take : data 

    main : 
        EXPANSION_HUNTER( data, params.genome )
}

/*
////////////////////////////////////////////////////////////////////
Workflow alternative entry point snpeff
////////////////////////////////////////////////////////////////////
*/

workflow snpeff_workflow {

    params.date = new Date().format( 'yyyyMMdd' )
    vcfs = Channel.fromPath( "${params.vcfs}/*.vcf" )
           
    snpeff( vcfs )

}

workflow post_variant_calling_gatk_workflow { 

    vcfs = Channel.fromPath( "${params.vcfs}/*.vcf.gz" )
    //vcfs = all_vcfs.take( 3 )

    // get the contigs 
    get_contigs( params.genome )

    // prepare genome
    gatk_prepare_genome( params.genome )

    // genotype cohort
    gatk_genotype_cohort( vcfs, gatk_prepare_genome.out.refs, get_contigs.out.contigs )
    genotype_vcfs = gatk_genotype_cohort.out.vcf_tbi

    // filters
    gatk_filters( genotype_vcfs, gatk_prepare_genome.out.refs, get_contigs.out.contig_file )

    // snpeff
    snpeff_input = gatk_genotype_cohort.out.vcfs.mix(gatk_filters.out.vcfs).flatten()
    snpeff( snpeff_input )

}


/*
////////////////////////////////////////////////////////////////////
Workflow
////////////////////////////////////////////////////////////////////
*/
workflow {

    if ( params.input_type == 'bam' ) {
        
        if (params.download_bams) {
            bam_list = Channel
                .fromPath( params.design )
                .splitCsv(header:false, sep : ',', skip : 1)
                .map { 
                    row -> [ row[0], row[1]] 
                    }
            
            download_bam( bam_list )
            raw_dat = download_bam.out.bams

        } else {
            raw_dat = Channel
                .fromPath( "${params.design}/*.bam" )
                .map { bam -> [ bam.SimpleName, bam] }

        }
        
        index_bam( raw_dat )
        bams = index_bam.out.bams

    } else {

        read_data_fq( params.design )

        raw_dat = read_data_fq.out.reads

        if (params.trim_adapters) { 
            trim( raw_dat ) 
            data = trim.out.fqs
        } else { 
            data = raw_dat         
        }

        bwa_align( data )

        bams = bwa_align.out.bams
    }

    if ( params.filter_bam ) { 
        filter_bam( bams ) 
        processed_bam = filter_bam.out.filter_bam_bai
    } else { 
        processed_bam = bams 
    }

    if ( params.call_variants ) {

        if ( params.variant_caller == 'DeepVariant' ) { 

            deepvariant( processed_bam )
            vcfs = deepvariant.out.vcfs
            snpeff_input = vcfs
        
        } else {

            gatk( processed_bam )

            if (params.genotype_cohort) {

                gatk_genotype_cohort( gatk.out.vcfs, gatk.out.refs, gatk.out.contigs )
                vcfs = gatk_genotype_cohort.out.vcf_tbi

            } else { 
                vcfs = gatk.out.vcfs
            }

            if (params.gatk_filters) {

                gatk_filters( vcfs, gatk.out.refs, gatk.out.contig_file )
                filt_vcf = gatk_filters.out.vcfs

            } else { 

                filt_vcf = vcfs
            
            }

            snpeff_input = gatk_genotype_cohort.out.vcfs.mix(gatk_filters.out.vcfs).flatten()
        }

        if ( params.run_snpeff ) {
            snpeff( snpeff_input )
        }
    }

    if (params.run_expansion_hunter) {

        expansion_hunter( processed_bam )   
    
    }

}