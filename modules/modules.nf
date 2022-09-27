process CHIP_SAMPLES_SHEET {

    label 'local'

    publishDir "$params.results/samples", mode : 'copy', pattern : '*'

    input : 
        val(samples)

    output : 
        path("controls.csv"), emit : control_ch
        path("fastq.csv"), emit : fq_ch
        path("replicates.csv"), emit : replicates_ch

    shell : 
    """
    #!/bin/bash

    source activate rnaseq

    which bwa

    python3 ${params.bin}/process_design_input.py -input ${samples}

    """
}

process TRIM_GALORE {

    tag "${sampleID}_TrimGalore"

    label 'medium'
    
    publishDir "$params.results/trimmed_fastq", mode : 'copy', pattern : '*.fq.gz'
    publishDir "$params.results/trimmed_fastq", mode : 'copy', pattern : '_trimming_report.txt'

    input : 
        tuple val(sampleID), val(fastq)

    output :
        tuple val(sampleID), path("${sampleID}*.fq.gz"), emit : fq_ch

    script :
    fastq_command = params.single_end ? "${fastq}" : "${fastq[0]} ${fastq[1]}"
    paired = params.single_end ? "" : "--paired"
    """
    #!/bin/bash
    
    source activate rnaseq

    trim_galore -j ${task.cpus} ${paired} --gzip --basename ${sampleID} ${fastq_command}

    """

}

process FASTQC {

    tag "${sampleID}_FASTQC"

    label 'medium'

    publishDir "$params.results/fastqc", mode : 'copy', pattern : '*_fastqc.{zip,html}'

    input : 
        tuple val(sampleID), val(fastq)
    
    output :
        path("*_fastqc.{zip,html}")

    script :
    fastq_command = params.single_end ? "${fastq}" : "${fastq[0]} ${fastq[1]}"
    """
    #!/bin/bash

    source activate rnaseq

    fastqc ${fastq_command} -o \$PWD

    """
}

process STAR_INDEX {

    tag "${genome}_STAR_index"

    label 'high'

    publishDir "$params.star_index", mode : 'copy'
    
    input : 
        val genome
        val gtf
    
    output : 
        path("*")
        val("${params.star_index}"), emit : star_idx_ch

    script : 
    """
    #!/bin/bash

    source activate rnaseq
    
    STAR --runMode genomeGenerate \\
         --genomeDir \$PWD \\
         --genomeFastaFiles ${genome} \\
         --sjdbGTFfile ${gtf} \\
         --runThreadN ${task.cpus}
    """
}

process STAR_ALIGN {

    tag "${sampleID}_STAR_align"

    label 'high'

    publishDir "$params.results/star_alignments/bam", mode : 'copy', pattern : '*.bam'
    publishDir "$params.results/star_alignments/bam", mode : 'copy', pattern : '*.bai'
    publishDir "$params.results/star_alignments/logs", mode : 'copy', pattern : '*_Log.final.out'

    input : 
        val star_idx
        val gtf
        tuple val(sampleID), val(fastq)
    
    output :
        tuple val(sampleID), path("${sampleID}_Aligned.sortedByCoord.out.bam"), path("${sampleID}_Aligned.sortedByCoord.out.bam.bai"), emit : star_bams_ch
        path("*_Log.final.out")

    script : 
    fastq_command = params.single_end ? "${fastq}" : "${fastq[0]} ${fastq[1]}"
    """
    #!/bin/bash
    
    source activate rnaseq

    STAR --genomeDir ${star_idx} \\
         --runThreadN ${task.cpus} \\
         --outFilterMultimapNmax 999 \\
         --outFilterType BySJout \\
         --readFilesCommand zcat \\
         --alignIntronMax 10000 \\
         --outFileNamePrefix ${sampleID}_ \\
         --outFilterMismatchNoverReadLmax 0.04 \\
         --readFilesIn ${fastq_command} \\
         --outSAMtype BAM SortedByCoordinate \\
         --outFilterIntronMotifs RemoveNoncanonical
    
    samtools index -@ 12 ${sampleID}_Aligned.sortedByCoord.out.bam
    """
}

process BWA_INDEX {

    tag "${genome}_BWA_INDEX"

    label 'high'

    publishDir "$params.bwa_index", mode : 'copy'
    
    input : 
        val genome

    output : 
        path("*")
        val("${params.bwa_index}"), emit : bwa_index_ch

    script : 
    """
    #!/bin/bash

    source activate rnaseq

    prefix=\$(basename ${genome} .fa)
    bwa index -p \$prefix ${genome}

    """
}


process BWA_MEM {

    tag "${sampleID}_BWA_MEM"

    label 'high'

    publishDir "$params.results/bwa_alignments/bam", mode : 'copy', pattern : '*sorted.bam'
    publishDir "$params.results/bwa_alignments/bam", mode : 'copy', pattern : '*.bai'
    publishDir "$params.results/bwa_alignments/logs", mode : 'copy', pattern : '*.summary.txt'

    input : 
        val bwa_idx
        tuple val(sampleID), val(fastq)

    output : 
        tuple val(sampleID), path("*sorted.bam"), emit : bwa_bam_ch
        path("*.summary.txt")
    
    script :
    fastq_command = params.single_end ? "${fastq}" : "${fastq[0]} ${fastq[1]}" 
    """
    #!/bin/bash

    source activate rnaseq

    bwa mem ${bwa_idx} \\
        -t ${task.cpus} \\
        ${fastq_command} | samtools view -@ ${task.cpus} -b -h -f 0x2 -F 0x0100 -O BAM -o ${sampleID}.bam -

    samtools sort  -@ ${task.cpus} ${sampleID}.bam -T ${sampleID} -o ${sampleID}.sorted.bam

    samtools index -@ ${task.cpus} ${sampleID}.sorted.bam

    samtools flagstat ${sampleID}.sorted.bam > ${sampleID}.sorted.bam.summary.txt
    """
}

process MERGE_BAM {

    tag "${condition}_merge_bams"

    label 'high'

    publishDir "$params.results/bwa_alignments/merged_bam/bam", mode : 'copy', pattern : '*.bam'
    publishDir "$params.results/bwa_alignments/merged_bam/bam", mode : 'copy', pattern : '*.bai'
    publishDir "$params.results/bwa_alignments/merged_bam/stats", mode : 'copy', pattern : '*stats'
    publishDir "$params.results/bwa_alignments/merged_bam/stats", mode : 'copy', pattern : '*txt'

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

process FILTER_MERGE_BAM {

    tag "${condition}_filter_merge_bam"

    label 'high'

    publishDir "$params.results/bwa_alignments/filtered_bam/bam", mode : 'copy', pattern : '*.bam'
    publishDir "$params.results/bwa_alignments/filtered_bam/bam", mode : 'copy', pattern : '*.bai'
    publishDir "$params.results/bwa_alignments/filtered_bam/stats", mode : 'copy', pattern : '*stats'

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
    """
    #!/bin/bash

    source activate rnaseq

    name=\$(basename $bam .bam)

    samtools view \\
        $filter_params \\
        $dup_params \\
        $multimap_params \\
        $blacklist_params \\
        -b ${bam} > \$name.filtered.bam

    samtools sort  -@ ${task.cpus} \$name.filtered.bam -T \$name -o \$name.filtered.sorted.bam
    samtools index \$name.filtered.sorted.bam
    samtools flagstat \$name.filtered.sorted.bam > \$name.filtered.sorted.flagstat
    samtools idxstats \$name.filtered.sorted.bam > \$name.filtered.sorted.idxstats
    samtools stats \$name.filtered.sorted.bam > \$name.filtered.sorted.stats
    """   
}

process PICARD_METRICS {

    tag "${sampleID}_PICARD_METRICS"

    label 'high'

    publishDir "$params.results/picard_metrics", mode : 'copy', pattern : '*'

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

process SALMON_INDEX {

    tag "${transcripts}_SALMON_index"

    label 'medium'

    publishDir "$params.salmon_index", mode : 'copy', pattern : "*"

    input : 
        val transcripts
    
    output : 
        path("*")
        val("${params.salmon_index}"), emit : salmon_idx_ch
    
    script : 
    """
    #!/bin/bash

    source activate rnaseq

    salmon index --threads ${task.cpus} -t ${transcripts} --index \$PWD

    """

}

process SALMON {

    tag "${sampleID}_SALMON_map"

    label 'high'

    publishDir "$params.results/counts/salmon", mode : 'copy', pattern : "*.counts"

    input : 
        val salmon_idx
        tuple val(sampleID), val(fastq)

    output : 
        tuple val(sampleID), path("${sampleID}.salmon.counts"), emit : salmon_quant_ch
    
    script : 
    fastq_command = params.single_end ? "-r ${fastq}" : "-1 ${fastq[0]} -2 ${fastq[1]}"
    """
    #!/bin/bash

    source activate rnaseq

    salmon quant \\
        -i ${salmon_idx} \\
        --threads ${task.cpus} \\
        --libType SR \\
        ${fastq_command} \\
        -o ${sampleID}
    
    cat ${sampleID}/quant.sf | cut -f1,5 > ${sampleID}.salmon.counts
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
    """
    #!/bin/bash

    source activate DeepTools

    name=\$(basename ${bam} .bam)

    bamCoverage -b ${bam} \\
        -o \$name.bw \\
        -p ${task.cpus} \\
        -bs 5 \\
        --smoothLength 10 \\
        --normalizeUsing RPKM \\
        --exactScaling \\
        --filterRNAstrand forward

    """
}

process MACS2 {

    tag "${ip}_vs_${control}_MACS2"

    label 'medium'

    publishDir "$params.results/MACS2/log", mode: 'copy', pattern : "*log"
    publishDir "$params.results/MACS2/peaks", mode: 'copy', pattern : "*narrowPeak"
    publishDir "$params.results/MACS2/summits", mode: 'copy', pattern : "*summits.bed"

    input : 
        tuple val(ip), path(ip_bam), val(control), path(control_bam)
    
    output : 
        path("*log")
        path("*Peak")
        path("*bed")

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
        -t ${ip_bam} \\
        -c ${control_bam} \\
        $broad \\
        -f $format \\
        -g $params.macs_gsize \\
        -n $ip \\
        $pileup \\
        $fdr \\
        $pvalue \\
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

process COUNTER {

    tag "${sampleID}_FeatureCounts"

    label 'medium'

    publishDir "$params.results/counts/feature_counts", mode: 'copy', pattern : "*.txt"
    publishDir "$params.results/counts/feature_counts", mode: 'copy', pattern : "*.summary"
    publishDir "$params.results/counts/feature_counts", mode: 'copy', pattern : "*.log"
    publishDir "$params.results/counts/feature_counts", mode: 'copy', pattern : "*.tsv"

    input :
        tuple val(sampleID), path(bam), path(bai)
        val gtf

    output :
        tuple val(sampleID), path("${sampleID}.counts.tsv"), emit : feature_counts_ch
        path("${sampleID}.feature_counts.txt")
        path("*.feature_counts.txt.summary")
        path("*.feature_counts.log")
    
    script :
    libtype = params.single_end ? "" : "-p"
    """
    #!/bin/bash

    source activate rnaseq

    featureCounts \\
        -t exon \\
        -g gene_id \\
        -s 2 \\
        -T ${task.cpus} \\
        -a ${gtf} \\
        -o ${sampleID}.feature_counts.txt \\
        ${libtype} \\
        ${bam} 2> ${sampleID}.feature_counts.log

    cat ${sampleID}.feature_counts.txt | grep -v "^#" | cut -f 1,7 | sed -e 's/.umi.dedup.bam//g' | sed -e 's/_Aligned.sortedByCoord.out.bam//g' | awk '(NR>1)' > ${sampleID}.counts.tsv

    cp ${sampleID}.feature_counts.txt.summary tmp.summary 
    rm ${sampleID}.feature_counts.txt.summary
    cat tmp.summary | sed -e 's/.umi.dedup.bam//g' | sed -e 's/_Aligned.sortedByCoord.out.bam//g' > ${sampleID}.feature_counts.txt.summary

    """

}

process MASTER_TABLE {

    label 'low'

    publishDir "$params.results/counts", mode : 'copy', pattern : "*.tsv"

    input :
        val project_name
        val counts
        val gene_ids
        
    output : 
        path("*.tsv"), emit : master_tables_ch
    
    script :
    """
    #!/bin/bash

    source activate rnaseq

    python3 ${params.bin}/make_master.py \\
        -f "${counts}" \\
        -o ${project_name} \\
        ${gene_ids}

    """
}

process RBIND_COUNTS {

    label 'low'

    publishDir "$params.results/counts", mode : 'copy', pattern : "*.tsv"

    input :
        tuple val(sampleID), path(feature_counts), path(salmon)

    output :
        path("*.tsv"), emit : counts_ch

    script : 
    """
    #!/bin/bash

    source activate rnaseq

    cat ${salmon} ${feature_counts} | grep -v "Name" | grep -v "Geneid" >> ${sampleID}.tsv
    """

}

process MULTI_QC {

    tag "${dir}_MULTIQC"

    label 'low'

    publishDir "$params.results", mode : 'copy'

    input : 
        val dir
    
    output : 
        path("multiqc_report.html")
    
    script : 
    """
    #!/bin/bash

    source activate Multiqc

    multiqc ${dir}
    """

}