process DOWNLOAD_BAM {

    label 'download_bam'

    tag "${file}"

    errorStrategy 'ignore'

    publishDir "$params.results/downloaded_bams", mode: 'copy', pattern : "*.bam"

    input : 
        tuple val(sample), val(download_link)

    output : 
        tuple val(sample), path("*.bam"), emit : bam_ch

    script :
    """
    #!/bin/bash

    wget -O ${sample}.bam ${download_link}
    """
}