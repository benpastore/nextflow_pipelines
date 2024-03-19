process BUILD_SNPEFF_DB {
    label "high"

    input:
        val genome
        val gtf

    output:
        path("*")
        path(genome), emit : genome

    script:
    """
    #!/bin/bash
    
    module load java/21.0.2

    db_name=\$(basename ${params.genome} .fa)
    db_dir=${params.snpEff_data}/data/\$db_name
    mkdir \$db_dir
    cp ${params.genome} \$db_dir/sequences.fa
    cp ${params.gtf} \$db_dir/sequences.fa
    echo "\$db_name : ${params.organism}" > temp.config
    cat ${params.snpEff_data}/snpEff.config temp.config 
    cd ${params.snpEff_data}
    
    java -jar ${params.snpEff_data}/snpEff.jar build -gft22 -v \$db_name

    """
}


process RUN_SNPEFF {
    label ""

    publishDir "$params.results/snp_Eff", mode : 'copy', pattern : '*'

    input:
        val vcf
        val genome

    output:
        path("*.vcf")


    script:
    """
    #!/bin/bash
    
    module load java/21.0.2
    
    db_name=\$(basename ${params.genome} .fa)
    vcf_base=\$(basename ${vcf} .vcf)

    zcat ${vcf} | cut -f1-8 > temp.vcf
    java -Xmx110g -jar ${params.snpEff_data}/snpEff.jar \
            -v \
            -no-downstream \
            -no-intergenic \
            -no-upstream \
            -onlyProtein \
            -ud 0 \
            \$db_name \
            temp.vcf > \$vcf_base.snpEff.vcf
            
    """

    
}