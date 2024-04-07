process PREPROCESS_SNPEFF {

    label 'local'

    input : 
        val genome
        val gtf
        val protein
        val cds
        val organism
    
    output :
        path("*")
        path("done.txt"), emit : done

    script:
    """
    #!/bin/bash

    db_name=\$(basename ${genome} .fa)
    db_dir=${params.snpEff_data}/data/\$db_name
    
    mkdir -p \$db_dir
    
    cp ${genome} \$db_dir/sequences.fa
    cp ${genome} ./sequences.fa

    if [[ "$organism" == *"elegans" ]]; then 
        cat ${protein} | sed -e 's/ wormpep=.*/.1/g' > \$db_dir/protein.fa
        cat ${cds} | sed -e 's/ gene=.*/.1/g' > \$db_dir/cds.fa
    else 
        python3 ${params.bin}/fix_human_cds_protein_annotation.py -f ${protein} -o \$db_dir/protein.fa
        python3 ${params.bin}/fix_human_cds_protein_annotation.py -f ${cds} -o \$db_dir/cds.fa
    fi

    cp ${gtf} \$db_dir/genes.gtf

    echo "\$db_name.genome : ${organism}" > temp.config
    
    cat ${params.snpEff_data}/snpEff.config temp.config > ${params.snpEff_data}/snpEff.config

    echo "Fin" > done.txt
    """
}

process BUILD_SNPEFF_DB {

    label "snpeff"

    input:
        val genome
        
    output:
        val("db_name.txt"), emit : db_name

    script:
    """
    #!/bin/bash

    db_name=\$(basename ${params.genome} .fa)    

    cd ${params.snpEff_data}
    
    java -jar ${params.snpEff_data}/snpEff.jar build -gtf22 -v \$db_name

    echo "\$db_name" > db_name.txt
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
    
    if [[ "${vcf}" == *".gz" ]]; then 
        vcfbase=\$(basename ${vcf} .vcf.gz)
        zcat ${vcf} | cut -f1-8 > tmp
    else 
        vcfbase=\$(basename ${vcf} .vcf)
        cat ${vcf} | cut -f1-8 > tmp
    fi

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