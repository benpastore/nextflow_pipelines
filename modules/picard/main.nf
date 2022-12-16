process MERGE_BAM {

    tag "${condition}_merge_bams"

    label 'extra_high'

    publishDir "$params.results/merged/bam", mode : 'copy', pattern : '*.bam'
    publishDir "$params.results/merged/bam", mode : 'copy', pattern : '*.bai'
    publishDir "$params.results/merged/stats", mode : 'copy', pattern : '*stats'
    publishDir "$params.results/merged/stats", mode : 'copy', pattern : '*txt'

    input : 
        tuple val(sampleIDs), path(bams), val(condition)

    output : 
        tuple val(condition), path("*.merged.sorted.bam"), path("*.merged.sorted.bam.bai"), emit : merged_bams_ch
        path("*stats")
        path("*.txt")
    
    script :
    avail_mem = task.memory ? "${task.memory.toGiga()}" : '5'
    if (bams.size() > 1) {
        """
        #!/bin/bash

        source activate rnaseq

        merged=${condition}.merged.sorted.bam
        #marked=${condition}.merged.sorted.marked.bam
        #metrics=${condition}.merged.sorted.marked.metrics.txt

        picard -Xmx${avail_mem}g MergeSamFiles \\
            ${'INPUT='+bams.join(' INPUT=')} \\
            OUTPUT=\$merged \\
            SORT_ORDER=coordinate \\
            VALIDATION_STRINGENCY=LENIENT \\
            TMP_DIR=tmp
        
        samtools index ${condition}.merged.sorted.bam

        #picard -Xmx${avail_mem}g MarkDuplicates \\
        #    INPUT=\$merged \\
        #    OUTPUT=\$marked \\
        #    ASSUME_SORTED=true \\
        #    REMOVE_DUPLICATES=false \\
        #    METRICS_FILE=\$metrics \\
        #    VALIDATION_STRINGENCY=LENIENT \\
        #    TMP_DIR=tmp
        
        #samtools index \$marked
        #samtools idxstats \$merged > \$merged.idxstats
        #samtools flagstat \$merged > \$merged.flagstats
        #samtools stats \$merged > \$merged.stats

        """
    } else {
        """
        #!/bin/bash

        source activate rnaseq

        file=${bams[0]}
        name=\$(basename $file .bam)
        merged=${name}.merged.sorted.bam

        cat $file > $merged

        samtools index $merged

        #marked=${condition}.merged.sorted.marked.bam    
        #metrics=${condition}.merged.sorted.marked.metrics.txt

        #picard -Xmx${avail_mem}g MarkDuplicates \\      
        #    INPUT=${bams[0]} \\
        #    OUTPUT=\$marked \\
        #    ASSUME_SORTED=true \\
        #    REMOVE_DUPLICATES=false \\
        #    METRICS_FILE=\$metrics \\
        #    VALIDATION_STRINGENCY=LENIENT \\
        #    TMP_DIR=tmp
        
        #samtools index \$marked
        #samtools idxstats \$marked > \$marked.idxstats
        #samtools flagstat \$marked > \$marked.flagstat
        #samtools stats \$marked > \$marked.stats
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