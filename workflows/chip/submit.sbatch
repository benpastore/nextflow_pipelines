#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PAS1473 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

nextflow /fs/ess/PCON0160/ben/pipelines/nextflow_pipelines/workflows/chip/main.nf --design /fs/ess/PCON0160/ben/projects/attf6_m01/miska_chip.csv --results /fs/ess/PCON0160/ben/projects/attf6_m01/miska_chip_results --single_end true -profile cluster --outprefix miska_chip

