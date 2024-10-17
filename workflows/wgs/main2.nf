#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.vcfs = '/fs/scratch/PCON0160/ben_cendr_genotype_cohort/vcfs'

vcfs = Channel
    .fromPath( "${params.vcfs}/*.vcf" )
    .concat(Channel.fromPath( "${params.vcfs}/*.vcf.gz" ))

vcfs.view()