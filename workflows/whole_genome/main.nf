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
/*
* This part is setting up for the default values of the input (we can always change these value in the command line 
* in the run.sh in the sam_cendr folder). However, we DO NOT want to change these params
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
Validate inputs (these are things that being specified in run.sh in projects/2024_CeNDR)
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


/*
////////////////////////////////////////////////////////////////////
Workflow
////////////////////////////////////////////////////////////////////
*/
date = new Date().format( 'yyyyMMdd' )
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

    if ( params.fastqc ){
        FASTQC( reads_ch )
    }


    // 3/14 maybe switch out bwa for STAR??
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
    // 

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
        // line 158 is important
        bam_to_bw_input_ch = BWA_MEM.out.bwa_bam_bai_ch
    }

    /*
     * BAM --> BW
     */
    BAM_TO_BW( bam_to_bw_input_ch )
    
    /* 
     * call variants 
     */
    if ( params.variant_caller == 'DeepVariant') {

        // DeepVariant
        variant_caller_input = bam_to_bw_input_ch
        DEEPVARIANT_CALL_VARIANTS( variant_caller_input, params.genome )

    } else {

        // general structure follows pipeline found at https://github.com/AndersenLab/wi-gatk
        // with some modifications to better fit my existing pipeline.

        /*
        * Get contigs
        */
        GET_CONTIGS( bam_to_bw_input_ch.first() )

        // get contigs: I, II, III, IV, V, X
        contigs = GET_CONTIGS.out.splitText { it.strip() }

        GATK_PREPARE_GENOME( params.genome )

        GATK_PROCESS_BAM( bam_to_bw_input_ch )

        /*
         * Variant caller combines individual chromosomes
         * to each vcf file and will call the variants for 
         * each chromosome independently. This is more efficient
         * as it allows calling variants for one vcf file as separate jobs (1 per chromosome)
         * the number of jobs will be: # vcf * # chromosomes
         */

        variant_caller_input = GATK_PROCESS_BAM
            .out
            .processed_bam_ch
            .combine( contigs )
            .combine( GATK_PREPARE_GENOME.out.processed_genome )

        //println "Contigs:"
        contigs.view()
        //println "Prepared Genome:"
        //GATK_PREPARE_GENOME.out.processed_genome.view()
        //println "Processed BAM:"
        //GATK_PROCESS_BAM.out.processed_bam_ch.view()
        //println "Variant caller input:"
       //variant_caller_input.view()

        // GATK
        GATK_CALL_VARIANTS( variant_caller_input )

        //GATK_CALL_VARIANTS.out.groupTuple().view()

        // concat vcfs, combine chromsome calls for each strain
        CONCAT_STRAIN_GVCFS( GATK_CALL_VARIANTS.out.groupTuple() )
        
        //CONCAT_STRAIN_GVCFS.out.concat_vcf_ch.view()
       
        // write sample map 
        MAKE_SAMPLE_MAP( CONCAT_STRAIN_GVCFS.out.concat_vcf_ch.collect().toList() )

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
            .combine( GATK_PREPARE_GENOME.out.processed_genome ) | \
            GATK_GENOTYPE_COHORT

        // combine vcf
        GATK_GENOTYPE_COHORT.out.genotype_cohort_gvcf_db_out.collect() | CONCATENATE_VCF
        
        CONCATENATE_VCF.out.vcf.combine(GATK_PREPARE_GENOME.out.processed_genome) | SOFT_FILTER
        SOFT_FILTER.out.soft_filter_vcf.combine(contigs) | HARD_FILTER

        // filters
        SOFT_FILTER
            .out
            .soft_vcf_stats.concat(
                HARD_FILTER.out.hard_vcf_stats
                )
                .collect() | MULTIQC
                
        MULTIQC
            .out
            .for_report
            .combine(SOFT_FILTER.out.soft_report)
            .combine(SOFT_FILTER.out.hard_vcf_stats)| HTML_REPORT

        snp_eff_input = CONCATENATE_VCF
            .out
            .vcf
            .mix(SOFT_FILTER.out.vcf, HARD_FILTER.out.vcf)
            .flatten()
        
        BUILD_SNPEFF_DB( params.genome, params.gtf )

        RUN_SNPEFF( snp_eff_input, BUILD_SNPEFF_DB.out.genome )

    }

    /*
    // SNP EFF *** add code for SnpEff https://pcingola.github.io/SnpEff/snpeff/introduction/ ***
    Step 1. building snp eff database
    https://pcingola.github.io/SnpEff/snpeff/build_db/#step-2-option-1-building-a-database-from-gtf-files

    1. Make directory for new genome
    2. put annotation files in directory (gtf file, genome reference)
    3. add information to genome configuration file

    cd /path/to/snpEff
    java -jar snpEff.jar build -gft22 -v c_elegans.WS283

    Step 2. Running SNPEff
    curr_dir=$PWD
    output="$curr_dir/WI.20220216.hard-filter_v2"
    sbatch="$curr_dir/sbatch"
    vcf="/fs/ess/PCON0160/CaeNDR/WI.20220216.hard-filter.vcf.gz"
    snpeff="/fs/ess/PCON0160/ben/pipelines/snpEff/snpEff.jar"
    db="c_elegans.WS283"

    java -Xmx110g -jar $snpeff \
        -v \
        -no-downstream \
        -no-intergenic \
        -no-upstream \
        -onlyProtein \
        -ud 0 \
        $db \
        $vcf \
        > $output/$name.ann.vcf


    */

    if ( params.run_expansion_hunter ) {
        // Repeat Expansion
        EXPANSION_HUNTER( variant_caller_input, params.genome )
    }

    // Indel caller
}
