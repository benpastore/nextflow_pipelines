process CI_SPLICE_AI {

    label "cispliceai"

    //publishDir "$params.results/CI_SPLICE", mode : 'copy', pattern : '*'

    input:
        tuple val(sampleID), val(vcf)
        val genome

    output:
        tuple val(sampleID), path("*.CISPLICEAI.vcf"), emit : ci_splice_ai_vcfs

    script :
    """
    #!/bin/bash

    module load cuda/12.3.0

    export TF_FORCE_UNIFIED_MEMORY='1'
    
    name=\$(basename ${vcf} .vcf)
    
    zcat ${vcf} > uncompressed.vcf
    
    singularity exec \\
        --nv \\
        --env LD_LIBRARY_PATH=\$LD_LIBRARY_PATH \\
        --env PATH=\$PATH \\
        --env TF_FORCE_GPU_ALLOW_GROWTH=true \\
        --env TF_GPU_ALLOCATOR=cuda_malloc_async \\
        --env TF_FORCE_UNIFIED_MEMORY='1' \\
        docker://benpasto/spliceai:latest \\
        cis-vcf --input uncompressed.vcf \\
        --output \${name}.CISPLICEAI.vcf \\
        -a grch38 ${genome}
    """
}

process CONCAT_CI_SPLICE_AI_FILES {

    label 'bcftools_low'

    publishDir "$params.results/CI_SPLICE", mode : 'copy', pattern : '*.CIspliceAI_combined.vcf'

    input : 
        tuple val(sampleID), val(vcfs)

    output : 
        tuple val(sampleID), val("*.CIspliceAI_combined.vcf"), emit : vcfs
    
    script : 
    """
    #!/bin/bash

    bcftools concat -o ${sampleID}.CIspliceAI_combined.vcf -O v ${vcfs.join(' ')}

    """
}