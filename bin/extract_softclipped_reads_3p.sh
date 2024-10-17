#!/bin/bash

module load python
source activate rnaseq_basic
cd /fs/ess/PCON0160/ben/projects/2024_july_dis3l2_human/hela_band/results_v4/20240717/star_alignments/bam

bam="/fs/ess/PCON0160/ben/projects/2024_july_dis3l2_human/hela_band/results_v4/20240717/star_alignments/bam/Hela-dis3l2-band_S17_R1_001.Aligned.sortedByCoord.out.bam"

#samtools view -h sorted.bam | awk '$6 ~ /S/' > soft_clipped_reads.sam
samtools view -h $bam |\
    awk '$6 ~ /S/' > soft_clipped_reads.sam

awk '{
  if ($6 ~ /[0-9]+S$/ || $6 ~ /^[0-9]+S/)
    print $0 
}' soft_clipped_reads.sam > soft_clipped_3p_reads.sam

