process BUILD_SNPEFF_DB {

    label "snpeff"

    input:
        val genome
        val gtf

    output:
        path("*")
        path("sequences.fa"), emit : genome

    script:
    """
    #!/bin/bash

    db_name=\\$(basename ${params.genome} .fa)
    db_dir=${params.snpEff_data}/data/\$db_name
    
    mkdir -p \$db_dir
    
    cp ${params.genome} \$db_dir/sequences.fa
    cp ${params.genome} ./sequences.fa
    cat ${params.protein} | sed -e 's/ wormpep=.*/.1/g' > \$db_dir/protein.fa
    cat ${params.cds} | sed -e 's/ gene=.*/.1/g' > \$db_dir/cds.fa
    cp ${params.gtf} \$db_dir/genes.gtf

    echo "\$db_name.genome : ${params.organism}" > temp.config
    
    cat ${params.snpEff_data}/snpEff.config temp.config > ${params.snpEff_data}/snpEff.config

    cd ${params.snpEff_data}
    
    java -jar ${params.snpEff_data}/snpEff.jar build -gtf22 -v \$db_name

    """
}


process RUN_SNPEFF {

    label "snpeff"

    publishDir "$params.results/SNPEFF", mode : 'copy', pattern : '*'

    input:
        val vcf
        val genome

    output:
        path("*.vcf")


    script:
    """
    #!/bin/bash
    
    db_name=\$(basename ${params.genome} .fa)

    vcfbase=\$(basename ${vcf} .vcf.gz)

    zcat ${vcf} | cut -f1-8 > tmp

    java -Xmx110g -jar ${params.snpEff_data}/snpEff.jar \
            -v \
            -no-downstream \
            -no-intergenic \
            -no-upstream \
            -onlyProtein \
            -ud 0 \
            \$db_name \
            tmp > \$vcfbase.snpEff.vcf
            
    """

    
}