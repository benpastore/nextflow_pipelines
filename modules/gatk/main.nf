
//https://github.com/AndersenLab/wi-gatk/blob/master/main.nf

process GATK_CALL_VARIANTS {

    tag "${condition}_filter_merge_bam"

    label 'GATK'

    publishDir "$params.results/gatk", mode : 'copy', pattern : '*'

    input : 
        tuple val(condition), val(bam), val(bai)
        val genome

    output : 
        tuple val(condition), val("*.vcf"), emit : gatk_vcf_ch
        path("*")
    
    script:
    """
    #!/bin/bash

    gatk HaplotypeCaller --java-options "-Xmx${task.memory.toGiga()}g -Xms1g -XX:ConcGCThreads=${task.cpus}" \\
        -R ${genome} \\
        -I ${bam} \\
        -O ${condition}.gatk.vcf
    """   
}