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