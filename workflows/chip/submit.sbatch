#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PAS1473 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

nextflow /fs/ess/PCON0160/ben/pipelines/nextflow_pipelines/workflows/chip/main.nf --design /fs/ess/PCON0160/ben/projects/piRNA_biogenesis/samples_pe.csv --results /fs/ess/PCON0160/ben/projects/piRNA_biogenesis/chip_results_10102022 --single_end false -profile cluster --outprefix chip_analysis_wen -resume

