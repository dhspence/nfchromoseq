process GET_COVERAGE {
    tag "$meta.id"
    label 'process_higher'
    label 'final_output'
    container "ghcr.io/dhslab/docker-python3:231224"

    input:
    tuple val(meta), path(align_output)
    tuple val(chromoseq_inputs), path("*")
    val(chromoseq_parameters)

    output:
    tuple val(meta), path("${meta.id}.coverage_report.tsv"), emit: report
    path("versions.yml")                                   , emit: versions

    script:
    """
    cat ${chromoseq_inputs.gene_bed} ${chromoseq_inputs.sv_bed} | gunzip -c | sort -k 1,1V -k 2,2n > regions.bed && \\
    get_dragen_coverage.py -r ${params.fasta} -c "${chromoseq_parameters.coverage_qc_levels}" -o ${meta.id}.coverage_report.tsv regions.bed ${meta.id}_tumor.cram

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_dragen_coverage.py: \$(get_dragen_coverage.py --version)
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}.coverage_report.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        get_dragen_coverage.py: \$(get_dragen_coverage.py --version)
    END_VERSIONS
    """
}
