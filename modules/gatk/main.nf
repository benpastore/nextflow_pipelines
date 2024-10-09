//https://github.com/AndersenLab/wi-gatk/blob/master/main.nf

/*
process GET_CONTIGS {

    label 'low'
    
    input:
        tuple val(conditions), path(bam), path(bai)

    output:
        path("contigs.txt"), emit : contigs

    script:
    """
    #!/bin/bash

    source activate rnaseq
    
    samtools idxstats ${bam} | cut -f 1 | grep -v '*' > contigs.txt
    """
}
*/


process GATK_PREPARE_GENOME {

    tag "prepare_genome"

    label 'GATK'

    input : 
        val genome
        val targets

    output : 
        tuple path("genome.fa"), path("genome.fa.fai"), path("genome.dict"), emit : processed_genome
        path("contigs.txt"), emit : contigs

    script:
    """
    #!/bin/bash

    cp ${genome} ./genome.fa

    gatk CreateSequenceDictionary -R ./genome.fa
    
    samtools faidx ./genome.fa

    python3 ${params.bin}/get_contigs.py -fai ./genome.fa.fai -target_contigs "${targets}"

    """
}

process IMPORT_GENOME_DB_CONTIGS {

    tag "prepare_genome"

    label 'GATK_xs'

    input : 
        val genome
        val targets

    output : 
        path("contigs.txt"), emit : contigs

    script:
    """
    #!/bin/bash

    cp ${genome} ./genome.fa
    
    samtools faidx ./genome.fa

    python3 ${params.bin}/get_contigs.py -fai ./genome.fa.fai -target_contigs "${targets}" -splits 3

    """
}

process GET_CONTIGS {

    label 'low'

    input : 
        val(genome)
    
    output: 
        path("contigs.txt"), emit : contigs
    
    script: 
    """
    #!/bin/bash

    source activate rnaseq

    cp ${genome} ./genome.fa

    samtools faidx ./genome.fa

    python3 ${params.bin}/get_contigs.py -fai ./genome.fa.fai

    """
}

process GATK_PROCESS_BAM { 

    tag "${condition}_process_bam"

    label 'GATK'

    input : 
        tuple val(condition), val(bam), val(bai)
    
    output : 
        tuple val(condition), path("*.dups.grouped.bam"), path("*.dups.grouped.bam.bai"), emit : processed_bam_ch

    script : 
    """
    #!/bin/bash

    name=\$(basename ${bam} .bam)

    gatk MarkDuplicates -I ${bam} -O \$name.dups.bam -M \$name.dups.txt

    gatk AddOrReplaceReadGroups -I \$name.dups.bam -O \$name.dups.grouped.bam -RGID ${condition} -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 20
    
    samtools index \$name.dups.grouped.bam

    """


}

process GATK_CALL_VARIANTS {

    errorStrategy 'retry'
    maxRetries 2
    time { 2.hour * task.attempt } 

    tag "${condition}_filter_merge_bam"

    label 'GATK_low'

    //publishDir "$params.results/gatk", mode : 'copy', pattern : '*'

    input : 
        tuple val(condition), val(bam), val(bai), val(region), val(genome_fa), val(genome_index), val(genome_dict)

    output : 
        tuple val(condition), path("*.sorted.g.vcf.gz"), emit : gatk_vcf_ch
        //path("*.g.vcf"), emit : gatk_vcf_only_ch
        //path("*")       
    
    script:
    // -Xms1g -XX:ConcGCThreads=${task.cpus}
    """
    #!/bin/bash

    gatk HaplotypeCaller --java-options "-Xmx${task.memory.toGiga()}g " \\
        --emit-ref-confidence GVCF \\
        --annotation DepthPerAlleleBySample \\
        --annotation Coverage \\
        --annotation GenotypeSummaries \\
        --annotation TandemRepeat \\
        --annotation StrandBiasBySample \\
        --annotation ChromosomeCounts \\
        --annotation ReadPosRankSumTest \\
        --annotation AS_ReadPosRankSumTest \\
        --annotation AS_QualByDepth \\
        --annotation AS_StrandOddsRatio \\
        --annotation AS_MappingQualityRankSumTest \\
        --annotation DepthPerSampleHC \\
        --annotation-group StandardAnnotation \\
        --annotation-group AS_StandardAnnotation \\
        --annotation-group StandardHCAnnotation \\
        --do-not-run-physical-phasing \\
        -R ${genome_fa} \\
        -I ${bam} \\
        -L ${region} \\
        -O ${condition}.g.vcf

    bcftools sort ${condition}.g.vcf -o ${condition}.sorted.g.vcf
    bgzip ${condition}.sorted.g.vcf
    tabix ${condition}.sorted.g.vcf.gz

    """   
}

process  CONCAT_STRAIN_GVCFS {

    label 'bcftools'
    
    tag "${condition}"

    publishDir "$params.results/gatk", mode : 'copy', pattern : '*vcf*'

    input:
        tuple val(condition), val(vcfs)
    
    output:
        tuple path("${condition}.g.vcf.gz"), path("${condition}.g.vcf.gz.tbi")
        //path("${condition}.g.vcf.gz"), emit : concat_vcf_ch
        tuple val(condition), path("${condition}.g.vcf.gz"), emit : concat_vcf_ch

    script:
    """
    #!/bin/bash

    #ls *_cohort_pol.vcf.gz > contig_set.vcf.gz

    bcftools concat -a -Oz ${vcfs.join(' ')} > ${condition}.g.vcf.gz

    bcftools index --tbi ${condition}.g.vcf.gz
    
    """
}

process REMOVE_VCF_DUPS {

    label 'bcftools'

    tag "${condition}"

    publishDir  "$params.results/vcf_filtered", mode : 'copy', pattern : '*vcf*'

    input : 
        tuple val(condition), val(vcf)
    
    output : 
        path("${condition}.g.dups.vcf.gz"), emit : concat_vcf_ch
    
    script : 
    """
    #!/bin/bash

    bcftools norm -D -Oz ${vcf} > ${condition}.g.dups.vcf.gz

    bcftools index --tbi ${condition}.g.dups.vcf.gz

    """

}

process MAKE_SAMPLE_MAP {

    label 'local'

    tag 'make_sample_map'

    publishDir "$params.results/sample_map", mode : 'copy', pattern : 'sample_map.tsv'

    input : 
        val samples
    
    output : 
        path("sample_map.tsv"), emit : sample_map

    shell:
    """
    #!/usr/bin/env python3
    import os

    lines = ''
    with open("$samples", 'r') as f : 
        for line in f : 
            v = line.strip()
            lines += f"{os.path.basename(v).split('.')[0]}\\t{v}\\n"
    
    op = open('sample_map.tsv', 'w')
    op.write(lines)
    op.close()

    """
}

process IMPORT_GENOME_DB {

    label 'WIGATK'

    errorStrategy 'retry'
    time { 5.hour * task.attempt } 
    maxRetries 3
    
    tag "import_genome_db"

    input : 
        tuple val(sample_map), val(contig)

    output : 
        tuple val(contig), path("${contig}.db"), emit: genome_db_output

    script:
    """
    #!/bin/bash

    gatk --java-options "-Xmx${task.memory.toGiga()}g -XX:ConcGCThreads=${task.cpus}" \\
        GenomicsDBImport \\
        --genomicsdb-workspace-path ${contig}.db \\
        --batch-size 50 \\
        -L ${contig} \\
        --sample-name-map ${sample_map} \\
        --reader-threads ${task.cpus}

    """
}

process GATK_GENOTYPE_COHORT {

    label 'WIGATK_xl'

    tag 'genotype_cohort'

    input : 
        tuple val(contig), val(db), val(genome_fa), val(genome_index), val(genome_dict)

    output : 
        tuple path("${contig}_cohort_pol.vcf.gz"), path("${contig}_cohort_pol.vcf.gz.csi"), emit : genotype_cohort_gvcf_db_out

    script:
    """
    #!/bin/bash

    gatk --java-options "-Xmx${task.memory.toGiga()}g" GenotypeGVCFs \\
        -R ${genome_fa} \\
        -V gendb://${db} \\
        -G StandardAnnotation \\
        -G AS_StandardAnnotation \\
        -G StandardHCAnnotation \\
        -L ${contig} \\
        -O ${contig}_cohort.vcf

    if [ "${contig}" == "${params.mito_name}" ]
        then
            bcftools view -O b ${contig}_cohort.vcf > ${contig}_cohort.bcf
            bcftools index ${contig}_cohort.bcf
        else
            bcftools view -O z ${contig}_cohort.vcf > ${contig}_cohort.bcf
            #| ${params.bin}/het_polarization.nim
            bcftools index ${contig}_cohort.bcf
    fi

    bcftools view -O v --min-af 0.000001 --threads=${task.cpus-1} ${contig}_cohort.bcf | \\
    vcffixup - | \\
    bcftools view -O z > ${contig}_cohort_pol.vcf.gz
    bcftools index ${contig}_cohort_pol.vcf.gz
    """
}

process CONCATENATE_VCF {

    errorStrategy 'retry'
    maxRetries 3
    time { 5.hour * task.attempt } 
    
    label 'bcftools'

    publishDir "$params.results/vcfs", mode : 'copy', pattern : 'raw.sorted.vcf.gz'
    publishDir "$params.results/vcfs", mode : 'copy', pattern : 'raw.sorted.vcf.gz.tbi'


    input: 
      path("*")

    output:
        tuple path("raw.sorted.vcf.gz"), path("raw.sorted.vcf.gz.tbi"), emit: vcf_tbi
        path("raw.sorted.vcf.gz"), emit: 'vcf'

    """
    #!/bin/bash
        
    ls *_cohort_pol.vcf.gz > contig_set.vcf.gz

    #bcftools sort contig_set.vcf.gz -Oz -o contig_set.sorted.vcf.gz

    bcftools concat -Oz --file-list contig_set.vcf.gz > raw.vcf.gz
    
    bcftools sort raw.vcf.gz -Oz -o raw.sorted.vcf.gz

    bcftools index --tbi raw.sorted.vcf.gz

    """
}

process MULTIQC {

    publishDir "${params.results}/report", mode: 'copy'

    input:
        path("*")

    output:
        val("multiqc_data/*.json")
        path "multiqc.html", emit: "for_report"

    script:
    """
    #!/bin/bash
    
    multiqc -k json --filename multiqc.html .
    """

}

process SOFT_FILTER {

    label "WIGATK"

    publishDir "$params.results/vcfs", mode : 'copy', pattern : '*vcf*'
    
    input:
        tuple val(vcf), val(vcf_index), val(fa), val(fai), val(dict)

    output:
        tuple path("soft-filter.${params.date}.vcf.gz"), path("soft-filter.${params.date}.vcf.gz.csi"), emit: soft_filter_vcf
        path "soft-filter.${params.date}.vcf.gz.tbi"
        path "soft-filter.${params.date}.stats.txt", emit: 'soft_vcf_stats'
        tuple path("soft-filter.${params.date}.filter_stats.txt"), path("soft-filter.${params.date}.stats.txt"), emit: 'soft_report'
        path("*soft-filter*.vcf.gz"), emit: vcf
    
    script:
    """
    #!/bin/bash

    #function cleanup {
    #    rm out.vcf.gz
    #}
    #trap cleanup EXIT

    gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms1g" \\
        VariantFiltration \\
        -R ${fa} \\
        --variant ${vcf} \\
        --genotype-filter-expression "DP < ${params.min_depth}"    --genotype-filter-name "DP_min_depth" \\
        --filter-expression "QUAL < ${params.qual}"                --filter-name "QUAL_quality" \\
        --filter-expression "FS > ${params.fisherstrand}"          --filter-name "FS_fisherstrand" \\
        --filter-expression "QD < ${params.quality_by_depth}"      --filter-name "QD_quality_by_depth" \\
        --filter-expression "SOR > ${params.strand_odds_ratio}"    --filter-name "SOR_strand_odds_ratio" \\
        --genotype-filter-expression "isHet == 1"                  --genotype-filter-name "is_het" \\
        -O out.vcf
    
    bgzip out.vcf
    bcftools index --tbi out.vcf.gz
    
    # Apply high missing and high heterozygosity filters
    bcftools filter --threads ${task.cpus} --soft-filter='high_missing' --mode + --include 'F_MISSING  <= ${params.high_missing}' out.vcf.gz |\\
    bcftools filter --threads ${task.cpus} --soft-filter='high_heterozygosity' --mode + --include '( COUNT(GT="het") / N_SAMPLES ) <= ${params.high_heterozygosity}' -O z > soft-filter.${params.date}.vcf.gz
    bcftools index soft-filter.${params.date}.vcf.gz
    bcftools index --tbi soft-filter.${params.date}.vcf.gz
    bcftools stats --threads ${task.cpus} \\
                   -s- --verbose soft-filter.${params.date}.vcf.gz > soft-filter.${params.date}.stats.txt
    {
        echo -e 'QUAL\\tQD\\tSOR\\tFS\\tFILTER';
        bcftools query -f '%QUAL\t%INFO/QD\t%INFO/SOR\t%INFO/FS\t%FILTER\n' soft-filter.${params.date}.vcf.gz;
    }     > soft-filter.${params.date}.filter_stats.txt
    """
}

process HARD_FILTER {

    label "WIGATK"

    publishDir "$params.results/vcfs", mode : 'copy', pattern : '*vcf*'

    input:
        tuple val(vcf), val(vcf_index), val(contigs)

    output:
        tuple path("hard-filter.${params.date}.vcf.gz"), path("hard-filter.${params.date}.vcf.gz.csi"), emit: 'vcf'
        path "hard-filter.${params.date}.vcf.gz.tbi"
        path "hard-filter.${params.date}.stats.txt", emit: 'hard_vcf_stats'
        path("*hard-filter*.vcf.gz"), emit : hard_filter_vcf

    script:
    """
    #!/bin/bash

    #function cleanup {
    #    # cleanup files on completion
    #    rm I.vcf.gz II.vcf.gz III.vcf.gz IV.vcf.gz V.vcf.gz X.vcf.gz MtDNA.vcf.gz
    #}
    #trap cleanup EXIT

    # Generate hard-filtered VCF
    function generate_hard_filter {

        if [ \${1} == "${params.mito_name}" ]
            then
                bcftools view -m2 -M2 --trim-alt-alleles -O u --regions \${1} ${vcf} |\\
                bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT ~ "DP_min_depth"' |\\
                bcftools filter -O u --exclude 'FILTER != "PASS"' |\\
                bcftools view -O v --min-af 0.000001 --max-af 0.999999 |\\
                vcffixup - | \\
                bcftools view -O z --trim-alt-alleles > \${1}.vcf.gz

            else
                bcftools view -m2 -M2 --trim-alt-alleles -O u --regions \${1} ${vcf} |\\
                bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT ~ "DP_min_depth" | FORMAT/FT ~"is_het"' |\\
                bcftools filter -O u --exclude 'FILTER != "PASS"' |\\
                bcftools view -O v --min-af 0.000001 --max-af 0.999999 |\\
                vcffixup - | \\
                bcftools view -O z --trim-alt-alleles > \${1}.vcf.gz
        fi

    }

    export -f generate_hard_filter

    parallel --verbose generate_hard_filter :::: ${contigs}

    awk '{ print \$0 ".vcf.gz" }' ${contigs} > contig_set.tsv
    
    bcftools concat -O z --file-list contig_set.tsv > hard-filter.${params.date}.vcf.gz
    
    bcftools index hard-filter.${params.date}.vcf.gz
    
    bcftools index --tbi hard-filter.${params.date}.vcf.gz
    
    bcftools stats -s- --verbose hard-filter.${params.date}.vcf.gz > hard-filter.${params.date}.stats.txt

    """
}

process HTML_REPORT {

    container 'andersenlab/r_packages:latest'

    publishDir "$params.results/html_report", mode: 'copy'

    input:
        tuple path("multiqc.html"), path("soft_filter_filter"), path("soft_filter_stats"), path("hard_filter_stats")

    output:
        file("*.html")

    """

    cat "${workflow.projectDir}/bin/gatk_report.Rmd" | \\
        sed -e 's/RELEASE_DATE/${params.date}/g' |
        sed -e 's+variation/${params.date}.soft-filter.stats.txt+${soft_filter_stats}+g' |
        sed -e 's+variation/${params.date}.hard-filter.stats.txt+${hard_filter_stats}+g' |
        sed -e 's+variation/${params.date}.soft-filter.filter_stats.txt+${soft_filter_filter}+g' > gatk_report_${params.date}.Rmd

    Rscript -e "rmarkdown::render('gatk_report_${params.date}.Rmd')"
    """

}