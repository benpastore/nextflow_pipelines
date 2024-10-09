process FILTER_BAM {

    errorStrategy 'retry'
    maxRetries 3
    time { 5.hour * task.attempt } 
    
    tag "${condition}_filter_merge_bam"

    label 'high'

    publishDir "$params.results/filtered/bam", mode : 'copy', pattern : '*.bam'
    publishDir "$params.results/filtered/bam", mode : 'copy', pattern : '*.bai'
    publishDir "$params.results/filtered/log", mode : 'copy', pattern : '*stat'
    publishDir "$params.results/filtered/log", mode : 'copy', pattern : '*stat'

    input : 
        tuple val(condition), val(bam), val(bai)

    output : 
        tuple val(condition), path("*filtered.sorted.bam"), emit : filtered_bams_ch
        tuple val(condition), path("*filtered.sorted.bam"), path("*filtered.sorted.bam.bai"), emit : bam_to_bw_input
        path("*stats")
    
    script:
    filter_params = params.single_end ? '-F 0x004' : '-F 0x004 -F 0x0008 -f 0x001'
    dup_params = params.keep_dups ? '' : '-F 0x0400'
    multimap_params = params.keep_multi_map ? '' : '-q 20'
    blacklist_params = params.blacklist ? "-L $bed" : ''
    remove_orphan = params.remove_orphan && !params.single_end ? "-f 0x2" : ''
    select_primary = params.select_primary ? "-F 0x0100" : ''
    """
    #!/bin/bash

    source activate rnaseq

    name=\$(basename $bam .bam)

    samtools view \\
        $filter_params \\
        $dup_params \\
        $multimap_params \\
        $blacklist_params \\
        $remove_orphan \\
        $select_primary \\
        -h \\
        -b ${bam} > \$name.filtered.bam

    samtools sort  -@ ${task.cpus} \$name.filtered.bam -T \$name -o \$name.filtered.sorted.bam
    samtools index \$name.filtered.sorted.bam
    samtools flagstat \$name.filtered.sorted.bam > \$name.filtered.sorted.flagstat
    samtools idxstats \$name.filtered.sorted.bam > \$name.filtered.sorted.idxstats
    samtools stats \$name.filtered.sorted.bam > \$name.filtered.sorted.stats
    """   
}

process INDEX_BAM {

    tag "${condition}_index_bam"

    label 'high'

    errorStrategy 'ignore'

    publishDir "$params.results/bams", mode : 'copy', pattern : '*.bam*'

    input : 
        tuple val(condition), val(bam)

    output : 
        tuple val(condition), path("*.bam"), path("*bai")
    
    script :
    """
    #!/bin/bash

    source activate rnaseq

    cp ${bam} .
    name=\$(basename ${bam})
    samtools index \$name 
    
    """
}