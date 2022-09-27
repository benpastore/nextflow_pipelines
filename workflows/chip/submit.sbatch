#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PAS1473 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

nextflow /users/PAS1473/benpasto1/ben/pipelines/RNA/workflows/chip/main.nf -profile cluster --design /fs/ess/PCON0160/ben/projects/piRNA_biogenesis/samples_se.csv --results /fs/ess/PCON0160/ben/projects/piRNA_biogenesis/miska --single_end true

