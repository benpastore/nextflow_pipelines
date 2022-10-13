process BAM_COMPARE {

    tag "${ip}_vs_${control}_BAM_Compare"

    label 'high'

    publishDir "$params.results/BamCompare", mode: 'copy', pattern : "*.bw"

    input : 
        tuple val(ip), path(ip_bam), val(control), path(control_bam)

    output :
        path("*bw"), emit : bws_ch
    
    script :
    normalize = params.normalize_bw ? "--normalizeUsing ${params.normalize_bw} --scaleFactorsMethod None" : '--normalizeUsing RPKM --scaleFactorsMethod None'
    smooth = params.smooth_length ? "--smoothLength ${params.smooth_length}" : '--smoothLength 10' 
    exact = params.exact_scaling ? "--exactScaling" : ''
    filter_strand = params.filter_strand ? '--filterRNAstrand forward' : ''
    operation = params.bw_compare_op ? "--operation ${params.bw_compare_op}" : "--operation subtract"
    """
    #!/bin/bash

    source activate rnaseq

    samtools index -@ ${task.cpus} ${ip_bam}
    samtools index -@ ${task.cpus} ${control_bam}

    source activate DeepTools

    name=\$(basename ${ip_bam} .bam)

    bamCompare -b1 ${ip_bam} \\
        -b2 ${control_bam} \\
        -o ${ip}_vs_${control}_compare.bw \\
        ${normalize} \\
        ${smooth} \\
        ${exact} \\
        ${operation} \\
        ${filter_strand}
    """
}

process BAM_TO_BW {

    tag "${sampleID}_BAM_to_BW"

    label 'high'

    publishDir "$params.results/BigWig", mode: 'copy', pattern : "*.bw"

    input :
        tuple val(sampleID), path(bam), path(bai)

    output :
        tuple val(sampleID), path("*bw"), emit : bws_ch
    
    script :
    normalize = params.normalize_bw ? "--normalizeUsing ${params.normalize_bw}" : '--normalizeUsing RPKM'
    smooth = params.smooth_length ? "--smoothLength ${params.smooth_length}" : '--smoothLength 10' 
    exact = params.exact_scaling ? "--exactScaling" : ''
    filter_strand = params.filter_strand ? '--filterRNAstrand forward' : ''
    bin_size = params.bin_size ? "--binSize ${params.bin_size}" : '--binSize 10'
    """
    #!/bin/bash

    source activate DeepTools

    name=\$(basename ${bam} .bam)

    bamCoverage -b ${bam} \\
        -o \$name.bw \\
        -p ${task.cpus} \\
        ${normalize} \\
        ${smooth} \\
        ${exact} \\
        ${filter_strand} \\
        ${bin_size}
        
    """
}