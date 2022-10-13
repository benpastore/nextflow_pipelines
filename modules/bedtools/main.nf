process CHIP_INTERSECT {

    label 'high'

    publishDir "$params.results/BamIntersectPeaks", mode: 'copy', pattern : "*"

    input : 
        tuple val(sampleID), path(bam), path(peaks)

    output :
        path("*.peaks.bam"), emit : peaks_bam_ch
        path("*.peaks.bw"), emit : peaks_bws_ch
    
    script :
    normalize = params.normalize_bw ? "--normalizeUsing ${params.normalize_bw}" : '--normalizeUsing RPKM'
    smooth = params.smooth_length ? "--smoothLength ${params.smooth_length}" : '--smoothLength 10' 
    exact = params.exact_scaling ? "--exactScaling" : ''
    filter_strand = params.filter_strand ? '--filterRNAstrand forward' : ''
    bin_size = params.bin_size ? "--binSize ${params.bin_size}" : '--binSize 10'
    """
    #!/bin/bash

    source activate rnaseq

    name=\$(basename ${bam} .bam)

    bedtools intersect -abam ${bam} -b ${peaks} > \$name.peaks.bam

    samtools index -@ ${task.cpus} \$name.peaks.bam

    source activate DeepTools
    
    bamCoverage -b \$name.peaks.bam \\
        -o \$name.peaks.bw \\
        -p ${task.cpus} \\
        ${normalize} \\
        ${smooth} \\
        ${exact} \\
        ${filter_strand} \\
        ${bin_size}

    """
}