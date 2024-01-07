process FILTER_CNV {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-cleutils:2"

    input:
    tuple val(meta), path (seg_file), path(align_output)
    tuple val(chromoseq_inputs), path("*")
    val(chromoseq_parameters)

    output:
    tuple val(meta), path("${meta.id}.cnv_annotated.vcf.gz*", arity:'2'), emit: vcf
    path "versions.yml",    emit: versions

    script:
    """
    process_cnvs.py -g ${meta.sex} -s ${chromoseq_parameters.min_filter_cna_size} \\
    -f ${chromoseq_parameters.min_filter_cna_abundance} -c ${chromoseq_inputs.cytobands} \\
    ${seg_file} ${meta.id}.cnv.vcf.gz | \\
    bcftools sort - | bgzip -c > ${meta.id}.cnv_annotated.vcf.gz && \\
    tabix -p vcf ${meta.id}.cnv_annotated.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n 1 | cut -d ' ' -f 2)
        \$(process_cnvs.py --version)
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}.cnv_annotated.vcf.gz"
    touch "${meta.id}.cnv_annotated.vcf.gz.tbi"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n 1 | cut -d ' ' -f 2)
        \$(process_cnvs.py --version)
    END_VERSIONS
    """
}
