process MAKE_FASTQLIST {
    label 'process_tiny'
    container "ghcr.io/dhslab/docker-python3:231224"

    input:
    tuple val(id), path(inputs)

    output:
    tuple val(id), path('fastq_list.csv'), path('*.run_info.csv'), emit: fastqlist
    path 'versions.yml', emit: versions

    script:
    // if inputs is a directory, add -d to the arg, otherwise add -f
    def flag = inputs.isDirectory() ? "-d" : "-f"

    """
    prepare_fastqs.py -i ${id} ${flag} ${inputs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(prepare_fastqs.py -v)
    END_VERSIONS
    """

    stub:
    // if inputs is a directory, add -d to the arg, otherwise add -f
    def flag = inputs.isDirectory() ? "-d" : "-f"

    """
    prepare_fastqs.py -i ${id} ${input_args} ${inputs}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(prepare_fastqs.py -v)
    END_VERSIONS
    """

}
