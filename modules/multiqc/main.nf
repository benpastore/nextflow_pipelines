process MULTI_QC {

    tag "${dir}_MULTIQC"

    label 'low'

    publishDir "$params.results", mode : 'copy'

    input : 
        val dir
    
    output : 
        path("multiqc_report.html")
    
    script : 
    """
    #!/bin/bash

    source activate Multiqc

    multiqc ${dir}
    """

}