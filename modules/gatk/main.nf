
//https://github.com/AndersenLab/wi-gatk/blob/master/main.nf

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


process GATK_PREPARE_GENOME {

    tag "prepare_genome"

    label 'GATK'

    input : 
        val genome
    
    output : 
        tuple path("genome.fa"), path("genome.fa.fai"), path("genome.dict"), emit : processed_genome
    
    script:
    """
    #!/bin/bash

    cp ${genome} ./genome.fa

    gatk CreateSequenceDictionary -R ./genome.fa
    
    samtools faidx ./genome.fa
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

    gatk AddOrReplaceReadGroups -I \$name.dups.bam -O \$name.dups.grouped.bam -RGID 4 -RGLB lib1 -RGPL ILLUMINA -RGPU unit1 -RGSM 20
    
    samtools index \$name.dups.grouped.bam

    """


}

process GATK_CALL_VARIANTS {

    tag "${condition}_filter_merge_bam"

    label 'GATK'

    //publishDir "$params.results/gatk", mode : 'copy', pattern : '*'

    input : 
        tuple val(condition), val(bam), val(bai), val(region), val(genome_fa), val(genome_index), val(genome_dict)

    output : 
        tuple val(condition), path("*.vcf"), emit : gatk_vcf_ch
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
        -O ${condition}.${region}.g.vcf

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
        path("${condition}.g.vcf.gz"), emit : concat_vcf_ch

    script:
    """
    #!/bin/bash

    #awk '{ print \$0 ".g.vcf.gz" }' > contig_set.tsv

    bcftools concat  -O z ${vcfs.join(' ')} > ${condition}.g.vcf.gz

    bcftools index --tbi ${condition}.g.vcf.gz
    
    """
}

process MAKE_SAMPLE_MAP {

    label 'local'

    tag 'make_sample_map'

    publishDir "$params.results/sample_map", mode : 'copy', pattern : 'sample_map.tsv'

    input : 
        val vcfs
    
    output : 
        path("sample_map.tsv"), emit : sample_map

    shell:
    """
    #!/usr/bin/env python3
    import os

    vcfs=str('${vcfs[0].join(' ')}').split(" ")
    print(vcfs)

    lines = ''
    for v in vcfs : 
        lines += f"{os.path.basename(v).split('.')[0]}\\t{v}\\n"
    
    op = open('sample_map.tsv', 'w')
    op.write(lines)
    op.close()

    """
}

process IMPORT_GENOME_DB {

    label 'WIGATK'
    
    tag "import_genome_db"

    input : 
        tuple val(sample_map), val(contig)

    output : 
        tuple val(contig), path("${contig}.db"), emit: genome_db_output

    script:
    """
    #!/bin/bash

    gatk \\
        GenomicsDBImport \\
        --genomicsdb-workspace-path ${contig}.db \\
        --batch-size 16 \\
        -L ${contig} \\
        --sample-name-map ${sample_map} \\
        --reader-threads ${task.cpus}

    """
}

process GATK_GENOTYPE_COHORT {

    label 'WIGATK'

    tag 'genotype_cohort'

    input : 
        tuple val(contig), val(db), val(genome_fa), val(genome_index), val(genome_dict)

    output : 
        tuple path("${contig}_cohort_pol.vcf.gz"), path("${contig}_cohort_pol.vcf.gz.csi"), emit : genotype_cohort_gvcf_db_out

    script:
    """
    #!/bin/bash

    gatk GenotypeGVCFs \\
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

    label 'bcftools'

    publishDir "$params.results/vcfs", mode : 'copy', pattern : '*vcf*'

    input: 
      path("*")

    output:
        tuple path("raw.vcf.gz"), path("raw.vcf.gz.tbi"), emit: 'vcf'

    """
    #!/bin/bash
        
    ls *_cohort_pol.vcf.gz > contig_set.tsv

    bcftools concat  -O z --file-list contig_set.tsv > raw.vcf.gz
    
    bcftools index --tbi raw.vcf.gz

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
        tuple path("soft-filter.${date}.vcf.gz"), path("soft-filter.${date}.vcf.gz.csi"), emit: soft_filter_vcf
        path "soft-filter.${date}.vcf.gz.tbi"
        path "soft-filter.${date}.stats.txt", emit: 'soft_vcf_stats'
        tuple path("soft-filter.${date}.filter_stats.txt"), path("soft-filter.${date}.stats.txt"), emit: 'soft_report'
        path("*.vcf.gz"), emit: vcf
    
    script:
    """
    #!/bin/bash

    function cleanup {
        rm out.vcf.gz
    }
    trap cleanup EXIT
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
    bcftools filter --threads ${task.cpus} --soft-filter='high_heterozygosity' --mode + --include '( COUNT(GT="het") / N_SAMPLES ) <= ${params.high_heterozygosity}' -O z > soft-filter.${date}.vcf.gz
    bcftools index soft-filter.${date}.vcf.gz
    bcftools index --tbi soft-filter.${date}.vcf.gz
    bcftools stats --threads ${task.cpus} \\
                   -s- --verbose soft-filter.${date}.vcf.gz > soft-filter.${date}.stats.txt
    {
        echo -e 'QUAL\\tQD\\tSOR\\tFS\\tFILTER';
        bcftools query -f '%QUAL\t%INFO/QD\t%INFO/SOR\t%INFO/FS\t%FILTER\n' soft-filter.${date}.vcf.gz;
    }     > soft-filter.${date}.filter_stats.txt
    """
}

process HARD_FILTER {
    label "WIGATK"

    publishDir "$params.results/vcfs", mode : 'copy', pattern : '*vcf*'

    input:
        tuple val(vcf), val(vcf_index), val(contigs)
    output:
        tuple path("hard-filter.${date}.vcf.gz"), path("hard-filter.${date}.vcf.gz.csi"), emit: 'vcf'
        path "hard-filter.${date}.vcf.gz.tbi"
        path "hard-filter.${date}.stats.txt", emit: 'hard_vcf_stats'
        path("*.vcf.gz"), emit: hard_filter_vcf

    script:
    """
    #!/bin/bash

    function cleanup {
            # cleanup files on completion
            rm I.vcf.gz II.vcf.gz III.vcf.gz IV.vcf.gz V.vcf.gz X.vcf.gz MtDNA.vcf.gz
        }
        trap cleanup EXIT

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
        bcftools concat  -O z --file-list contig_set.tsv > hard-filter.${date}.vcf.gz

        bcftools index hard-filter.${date}.vcf.gz
        bcftools index --tbi hard-filter.${date}.vcf.gz
        bcftools stats -s- --verbose hard-filter.${date}.vcf.gz > hard-filter.${date}.stats.txt
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
        sed -e 's/RELEASE_DATE/${date}/g' |
        sed -e 's+variation/WI.{vcf_date}.soft-filter.stats.txt+${soft_filter_stats}+g' |
        sed -e 's+variation/WI.{vcf_date}.hard-filter.stats.txt+${hard_filter_stats}+g' |
        sed -e 's+variation/WI.{vcf_date}.soft-filter.filter_stats.txt+${soft_filter_filter}+g' > gatk_report_${date}.Rmd
    Rscript -e "rmarkdown::render('gatk_report_${date}.Rmd')"
        
    """

}