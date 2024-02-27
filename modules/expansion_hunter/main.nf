process EXPANSION_HUNTER {

    tag "${sampleID}_EXPANSION_HUNTER"

    label 'expansion_hunter'
    
    publishDir "$params.results/expansion_hunter/${sampleID}", mode : 'copy', pattern : '*'

    input : 
        tuple val(sampleID), val(bam), val(bai)
        val genome

    output :
        path("*")


    script :
    """
    #!/bin/bash
    #singularity run docker://benpasto/atria:latest /run_atria -h

    /ExpansionHunterDenovo-v0.9.0-linux_x86_64/bin/ExpansionHunterDenovo profile --reads ${bam} --reference ${genome} --output-prefix ${sampleID}

    """

}
