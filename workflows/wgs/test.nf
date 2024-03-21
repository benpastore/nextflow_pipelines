#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.run = true

params.name = '/fs/ess/PCON0160/ben/genomes/c_elegans/WS279/c_elegans.PRJNA13758.WS279.genomic.fa'

workflow print_it {

    take : 
        data

    main :
        println params.run
        RETURN_NAME( params.name )

    emit : 
        dat = RETURN_NAME.out.output_ch

}

workflow {

    print_it( params.name )
    print_it.out.dat.view()

}

process RETURN_NAME {

    tag 'local'

    input : 
        val name
    
    output : 
        path("my_name.txt"), emit : output_ch
    
    script : 
    """
    #!/bin/bash

    which java
    java -h
    echo "Hello my name is ${name}" > my_name.txt
    """

}