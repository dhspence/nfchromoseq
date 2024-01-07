process FILTER_SV {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'

    container "ghcr.io/dhslab/docker-cleutils:2"
    
    input:
    tuple val(meta), path(vcf)
    tuple val(chromoseq_inputs), path("*")
    val(chromoseq_parameters)

    output:
    tuple val(meta), path("${meta.id}.sv_annotated.vcf.gz*", arity: '2'), emit: vcf
    path "versions.yml",    emit: versions

    script:
    """        
    filter_sv.py -l ${chromoseq_parameters.min_filter_sv_length} -L ${chromoseq_parameters.max_filter_sv_length} \\
    -g ${chromoseq_inputs.sv_targets} -m ${chromoseq_parameters.min_filter_sv_reads} \\
    -a ${chromoseq_parameters.min_filter_sv_abundance} -o ${meta.id}.sv_annotated.vcf.gz ${meta.id}.sv_vep.vcf.gz && \\
    tabix -p vcf ${meta.id}.sv_annotated.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(filter_sv.py --version)
    END_VERSIONS
    """

    stub: 
    """
    touch "${meta.id}.sv_annotated.vcf.gz"
    touch "${meta.id}.sv_annotated.vcf.gz.tbi"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(filter_sv.py --version)
    END_VERSIONS
    """
}