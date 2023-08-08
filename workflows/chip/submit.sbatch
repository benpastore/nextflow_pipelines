#!/bin/bash

#SBATCH --no-requeue
#SBATCH -c 1
#SBATCH --account=PCON0160 --nodes=1 --ntasks-per-node=1 --time=8:00:00 --mem=2gb

nextflow /fs/ess/PCON0160/ben/pipelines/nextflow_pipelines/workflows/chip/main.nf --design /fs/ess/PAS1473/ben_past_projects/Function_eggd_smRNA_RNA/samples_se_chip_hlh30_modencode.csv --results /fs/ess/PAS1473/ben_past_projects/Function_eggd_smRNA_RNA/hlh30_chip_modEncode --single_end true -profile cluster --outprefix analysis

