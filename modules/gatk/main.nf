
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
            -R ref.fa.gz \\
            -I ${bam} \\
            -L ${region} \\
            -O ${region}.gatk.vcf   

        bcftools view -O z ${region}.gatk.vcf > ${region}.gatk.vcf.gz
        
        rm ${region}.gatk.vcf
        
    """   
}

//https://github.com/AndersenLab/wi-gatk/blob/master/main.nf
process get_contigs {

    label 'sm'

    input:
        tuple strain, path(bam), path(bai)

    output:
        path("contigs.txt")

    """
        samtools idxstats ${bam} | cut -f 1 | grep -v '*' > contigs.txt
    """

}

process hard_filter {


    label 'lg'

    publishDir "${params.output}/variation", mode: 'copy'

    input:
        tuple path(vcf), path(vcf_index), path(contigs)

    output:
        tuple path("WI.${date}.hard-filter.vcf.gz"), path("WI.${date}.hard-filter.vcf.gz.csi"), emit: 'vcf'
        path "WI.${date}.hard-filter.vcf.gz.tbi"
        path "WI.${date}.hard-filter.stats.txt", emit: 'hard_vcf_stats'


    """
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
        bcftools concat  -O z --file-list contig_set.tsv > WI.${date}.hard-filter.vcf.gz

        bcftools index WI.${date}.hard-filter.vcf.gz
        bcftools index --tbi WI.${date}.hard-filter.vcf.gz
        bcftools stats -s- --verbose WI.${date}.hard-filter.vcf.gz > WI.${date}.hard-filter.stats.txt
    """
}

process soft_filter {

    publishDir "${params.output}/variation", mode: 'copy'

    label 'xl'

    input:
        tuple path(vcf), path(vcf_index), file("ref_file"), file("ref_index"), file("ref_dict"), file("ref_gzi")

    output:
        tuple path("WI.${date}.soft-filter.vcf.gz"), path("WI.${date}.soft-filter.vcf.gz.csi"), emit: soft_filter_vcf
        path "WI.${date}.soft-filter.vcf.gz.tbi"
        path "WI.${date}.soft-filter.stats.txt", emit: 'soft_vcf_stats'
        tuple path("WI.${date}.soft-filter.filter_stats.txt"), path("WI.${date}.soft-filter.stats.txt"), emit: 'soft_report'
    

    """
    cp ${ref_file} ./ref.fa.gz
    cp ${ref_index} ./ref.fa.gz.fai
    cp ${ref_dict} ./ref.dict
    cp ${ref_gzi} ./ref.fa.gz.gzi

        function cleanup {
            rm out.vcf.gz
        }
        trap cleanup EXIT

        gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms1g" \\
            VariantFiltration \\
            -R ref.fa.gz \\
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
        bcftools filter --threads ${task.cpus} --soft-filter='high_heterozygosity' --mode + --include '( COUNT(GT="het") / N_SAMPLES ) <= ${params.high_heterozygosity}' -O z > WI.${date}.soft-filter.vcf.gz

        bcftools index WI.${date}.soft-filter.vcf.gz
        bcftools index --tbi WI.${date}.soft-filter.vcf.gz
        bcftools stats --threads ${task.cpus} \\
                       -s- --verbose WI.${date}.soft-filter.vcf.gz > WI.${date}.soft-filter.stats.txt

        {
            echo -e 'QUAL\\tQD\\tSOR\\tFS\\tFILTER';
            bcftools query -f '%QUAL\t%INFO/QD\t%INFO/SOR\t%INFO/FS\t%FILTER\n' WI.${date}.soft-filter.vcf.gz;
        }     > WI.${date}.soft-filter.filter_stats.txt

    """
}

process genotype_cohort_gvcf_db {
    // Heterozygous polarization is also performed here.

    tag { "${contig}" }
    label 'xl'

    input:
        tuple val(contig), file("${contig}.db"), file("ref_file"), file("ref_index"), file("ref_dict"), file("ref_gzi")

    output:
        tuple file("${contig}_cohort_pol.vcf.gz"), file("${contig}_cohort_pol.vcf.gz.csi")


    /*
        het_polarization polarizes het-variants to REF or ALT (except for mitochondria)
    */

    """
    cp ${ref_file} ./ref.fa.gz
    cp ${ref_index} ./ref.fa.gz.fai
    cp ${ref_dict} ./ref.dict
    cp ${ref_gzi} ./ref.fa.gz.gzi

        gatk  --java-options "-Xmx${task.memory.toGiga()}g -Xms1g" \\
            GenotypeGVCFs \\
            -R ref.fa.gz \\
            -V gendb://${contig}.db \\
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
        
            bcftools view -O z ${contig}_cohort.vcf | \\
            het_polarization > ${contig}_cohort.bcf
            bcftools index ${contig}_cohort.bcf
        
        fi

        bcftools view -O v --min-af 0.000001 --threads=${task.cpus-1} ${contig}_cohort.bcf | \\
        vcffixup - | \\
        bcftools view -O z > ${contig}_cohort_pol.vcf.gz
        bcftools index ${contig}_cohort_pol.vcf.gz
    """
}

process concatenate_vcf {

    label 'lg'

    input: 
      path("*")

    output:
        tuple path("WI.raw.vcf.gz"), path("WI.raw.vcf.gz.tbi"), emit: 'vcf'

    """
        ls *_cohort_pol.vcf.gz > contig_set.tsv
        bcftools concat  -O z --file-list contig_set.tsv > WI.raw.vcf.gz
        bcftools index --tbi WI.raw.vcf.gz
    """
}

process soft_filter {

    publishDir "${params.output}/variation", mode: 'copy'

    label 'xl'

    input:
        tuple path(vcf), path(vcf_index), file("ref_file"), file("ref_index"), file("ref_dict"), file("ref_gzi")

    output:
        tuple path("WI.${date}.soft-filter.vcf.gz"), path("WI.${date}.soft-filter.vcf.gz.csi"), emit: soft_filter_vcf
        path "WI.${date}.soft-filter.vcf.gz.tbi"
        path "WI.${date}.soft-filter.stats.txt", emit: 'soft_vcf_stats'
        tuple path("WI.${date}.soft-filter.filter_stats.txt"), path("WI.${date}.soft-filter.stats.txt"), emit: 'soft_report'
    

    """
    cp ${ref_file} ./ref.fa.gz
    cp ${ref_index} ./ref.fa.gz.fai
    cp ${ref_dict} ./ref.dict
    cp ${ref_gzi} ./ref.fa.gz.gzi

        function cleanup {
            rm out.vcf.gz
        }
        trap cleanup EXIT

        gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms1g" \\
            VariantFiltration \\
            -R ref.fa.gz \\
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
        bcftools filter --threads ${task.cpus} --soft-filter='high_heterozygosity' --mode + --include '( COUNT(GT="het") / N_SAMPLES ) <= ${params.high_heterozygosity}' -O z > WI.${date}.soft-filter.vcf.gz

        bcftools index WI.${date}.soft-filter.vcf.gz
        bcftools index --tbi WI.${date}.soft-filter.vcf.gz
        bcftools stats --threads ${task.cpus} \\
                       -s- --verbose WI.${date}.soft-filter.vcf.gz > WI.${date}.soft-filter.stats.txt

        {
            echo -e 'QUAL\\tQD\\tSOR\\tFS\\tFILTER';
            bcftools query -f '%QUAL\t%INFO/QD\t%INFO/SOR\t%INFO/FS\t%FILTER\n' WI.${date}.soft-filter.vcf.gz;
        }     > WI.${date}.soft-filter.filter_stats.txt

    """
}