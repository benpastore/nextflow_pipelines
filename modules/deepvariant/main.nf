process DEEPVARIANT_CALL_VARIANTS {

    time { 5.hour * task.attempt } 
    errorStrategy 'retry'
    maxRetries 3 

    tag "${condition}_filter_merge_bam"

    label 'DeepVariant'

    publishDir "$params.results/deepvariant", mode : 'copy', pattern : '*'

    input : 
        tuple val(condition), val(bam), val(bai)
        val genome

    output : 
        tuple val(condition), val("*.vcf.gz"), emit : deepvariant_vcf_ch
        tuple val(condition), val("*.gvcf.gz"), emit : deepvariant_gvcf_ch
        path("*")
        path("*vcf.gz"), emit : vcf
    
    script:
    """
    #!/bin/bash

    export TF_FORCE_UNIFIED_MEMORY='1'

    singularity run \\
        docker://google/deepvariant:latest \\
        /opt/deepvariant/bin/run_deepvariant \\
        --model_type="WGS" \\
        --ref=${genome} \\
        --reads=${bam} \\
        --output_vcf=${condition}.dv.vcf \\
        --output_gvcf=${condition}.dv.gvcf \\
        --num_shards=${task.cpus} \\
        --logging_dir=\$PWD/${condition}_dv_logs

    gzip *vcf*
        
    """   
}