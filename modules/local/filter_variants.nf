process FILTER_VARIANTS {
    tag "$meta.id"
    label 'process_low'
    container "ghcr.io/dhslab/docker-cleutils:2"
     
    input:
    tuple val(meta), path(align_output)
    tuple val(chromoseq_inputs), path("*")
    val(chromoseq_parameters)

    output:
    tuple val(meta), path("${meta.id}.smallvariants.vcf.gz*", arity: '2'), emit: vcf
    path "versions.yml",    emit: versions

    script:
    """
    filter_variants.py -b ${chromoseq_inputs.gene_bed} -L ${chromoseq_parameters.max_smallvariant_sv_length} \\
    -m ${chromoseq_parameters.min_smallvariant_reads} -a ${chromoseq_parameters.min_smallvariant_vaf} \\
    ${meta.id}.hard-filtered.vcf.gz ${meta.id}.sv.vcf.gz | \\
    bcftools sort - | bgzip -c > ${meta.id}.smallvariants.vcf.gz && \\
    tabix -p vcf ${meta.id}.smallvariants.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_variants.py: \$(filter_variants.py --version)
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}.smallvariants.vcf.gz"
    touch "${meta.id}.smallvariants.vcf.gz.tbi"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        filter_variants.py: \$(filter_variants.py --version)
    END_VERSIONS
    """
}
