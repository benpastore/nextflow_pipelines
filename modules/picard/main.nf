process FILTER_MERGE_BAM {

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
    multimap_params = params.keep_multi_map ? '' : '-q 1'
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

process MERGE_BAM {

    tag "${condition}_merge_bams"

    label 'high'

    publishDir "$params.results/merged/bam", mode : 'copy', pattern : '*.bam'
    publishDir "$params.results/merged/bam", mode : 'copy', pattern : '*.bai'
    publishDir "$params.results/merged/stats", mode : 'copy', pattern : '*stats'
    publishDir "$params.results/merged/stats", mode : 'copy', pattern : '*txt'

    input : 
        tuple val(sampleIDs), path(bams), val(condition)

    output : 
        tuple val(condition), path("*.merged.sorted.marked.bam"), path("*.merged.sorted.marked.bam.bai"), emit : merged_bams_ch
        path("*stats")
        path("*.txt")
    
    script :
    avail_mem = task.memory ? "${task.memory.toGiga()}" : '5'
    if (bams.size() > 1) {
        """
        #!/bin/bash

        source activate rnaseq

        merged=${condition}.merged.sorted.bam
        marked=${condition}.merged.sorted.marked.bam
        metrics=${condition}.merged.sorted.marked.metrics.txt

        picard -Xmx${avail_mem}g MergeSamFiles \\
            ${'INPUT='+bams.join(' INPUT=')} \\
            OUTPUT=\$merged \\
            SORT_ORDER=coordinate \\
            VALIDATION_STRINGENCY=LENIENT \\
            TMP_DIR=tmp
        
        samtools index ${condition}.merged.sorted.bam

        picard -Xmx${avail_mem}g MarkDuplicates \\
            INPUT=\$merged \\
            OUTPUT=\$marked \\
            ASSUME_SORTED=true \\
            REMOVE_DUPLICATES=false \\
            METRICS_FILE=\$metrics \\
            VALIDATION_STRINGENCY=LENIENT \\
            TMP_DIR=tmp
        
        samtools index \$marked
        samtools idxstats \$marked > \$marked.idxstats
        samtools flagstat \$marked > \$marked.flagstats
        samtools stats \$marked > \$marked.stats

        """
    } else {
        """
        #!/bin/bash

        source activate rnaseq

        marked=${condition}.merged.sorted.marked.bam    
        metrics=${condition}.merged.sorted.marked.metrics.txt

        picard -Xmx${avail_mem}g MarkDuplicates \\      
            INPUT=${bams[0]} \\
            OUTPUT=\$marked \\
            ASSUME_SORTED=true \\
            REMOVE_DUPLICATES=false \\
            METRICS_FILE=\$metrics \\
            VALIDATION_STRINGENCY=LENIENT \\
            TMP_DIR=tmp
        
        samtools index \$marked
        samtools idxstats \$marked > \$marked.idxstats
        samtools flagstat \$marked > \$marked.flagstat
        samtools stats \$marked > \$marked.stats
        """
    }

}

process PICARD_METRICS {

    tag "${sampleID}_PICARD_METRICS"

    label 'high'

    publishDir "$params.results/picard/metrics", mode : 'copy', pattern : '*'

    input : 
        tuple val(sampleID), path(bam)
        path(fasta)
    
    output : 
        path("*metrics")
        path("*pdf")
    
    script:
    avail_mem = task.memory ? "${task.memory.toGiga()}" : '5'
    """
    #!/bin/bash

    source activate rnaseq

    name=\$(basename $bam .bam)

    picard -Xmx${avail_mem}g CollectMultipleMetrics \\
        INPUT=${bam} \\
        OUTPUT=\$name.CollectMultipleMetrics \\
        REFERENCE_SEQUENCE=$fasta \\
        VALIDATION_STRINGENCY=LENIENT \\
        TMP_DIR=tmp
    """
}