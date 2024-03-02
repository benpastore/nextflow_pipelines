
//https://github.com/AndersenLab/wi-gatk/blob/master/main.nf

process GATK_CALL_VARIANTS {

    tag "${condition}_filter_merge_bam"

    label 'GATK'

    publishDir "$params.results/gatk", mode : 'copy', pattern : '*'

    input : 
        tuple val(condition), val(bam), val(bai), val(contigs)
        val genome

    output : 
        tuple val(condition), val("*.vcf"), emit : gatk_vcf_ch
        path("*")
    
    script:
    // -Xms1g -XX:ConcGCThreads=${task.cpus}
    """
    #!/bin/bash

    cp ${genome} ./genome.fa
    gatk CreateSequenceDictionary -R ./genome.fa
    samtools faidx ./genome.fa

    name=\$(basename ${bam} .bam)

    gatk MarkDuplicates -I ${bam} -O \$name.mark_dups.bam -M \$name.mark_dups.txt

    gatk AddOrReplaceReadGroups -I \$name.mark_dups.bam -O \$name.mark_dups.groups.bam -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 20
    
    samtools index \$name.mark_dups.groups.bam


    gatk HaplotypeCaller --java-options "-Xmx${task.memory.toGiga()}g " \\
        -R ./genome.fa \\
        -I \$name.mark_dups.groups.bam \\
        -O ${condition}.gatk.vcf
    """   
}

process GET_CONTIGS {

    label 'low'
    
    input:
        tuple val(conditions), path(bam), path(bai)

    output:
        tuple val(conditions), path(bam), path(bai), path("contigs.txt"), emit : variant_caller_input

    script:
    """
    #!/bin/bash

    source activate rnaseq
    
    samtools idxstats ${bam} | cut -f 1 | grep -v '*' > contigs.txt
    """
}