process MACS2 {

    tag "${ip}_vs_${control}_MACS2"

    label 'medium'

    publishDir "$params.results/MACS2/peaks", mode: 'copy', pattern : "*narrowPeak"
    publishDir "$params.results/MACS2/peaks", mode: 'copy', pattern : "*Peak"
    publishDir "$params.results/MACS2/summits", mode: 'copy', pattern : "*summits*"

    input : 
        tuple val(control), val(ip), path(ip_bam), path(control_bam)
    
    output : 
        path("*Peak")
        path("*bed")
        tuple val(ip), path("${ip}_peaks.narrowPeak"), emit : macs2_peaks_ch

    script : 
    broad = params.narrow_peak ? '' : "--broad --broad-cutoff ${params.broad_cutoff}"
    format = params.single_end ? 'BAM' : 'BAMPE'
    pileup = params.save_macs_pileup ? '-B --SPMR' : ''
    fdr = params.macs_fdr ? "--qvalue ${params.macs_fdr}" : ''
    pvalue = params.macs_pvalue ? "--pvalue ${params.macs_pvalue}" : ''
    """
    #!/bin/bash

    source activate rnaseq

    macs2 callpeak \\
        -t ${ip_bam.join(' ')} \\
        -c ${control_bam.join(' ')} \\
        -f $format \\
        -g $params.macs_gsize \\
        -n $ip \\
        $pileup \\
        $fdr \\
        $pvalue \\
        --keep-dup all

    macs2 callpeak \\
        -t ${ip_bam.join(' ')} \\
        -c ${control_bam.join(' ')} \\
        -f $format \\
        -g $params.macs_gsize \\
        -n $ip \\
        $pileup \\
        $fdr \\
        $pvalue \\
        --broad --broad-cutoff 0.1 \\
        --keep-dup all

    # bed file output
    #1. chrom
    #2. start
    #3. end
    #4. name 
    #5. score
    #6. strand
    #7. signalValue
    #8. pValue
    #9. qValue
    #10. peak (point source called for this peak, 0 based offset from chromStart)
    """
}