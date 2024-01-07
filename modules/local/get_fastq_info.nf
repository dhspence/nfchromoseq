process GET_FASTQ_INFO {
    label 'process_tiny'
    container "ghcr.io/dhslab/docker-python3:231224"
 
    input:
    tuple val(id), path(read1), path(read2)

    output:
    tuple val(id), path('fastq_list.csv'), path('*.run_info.csv'), emit: fastqlist
    path 'versions.yml', emit: versions

    script: 
    """
    prepare_fastqs.py -i ${id} -1 ${read1} -2 ${read2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(prepare_fastqs.py -v)
    END_VERSIONS
    """

    stub:
    """
    prepare_fastqs.py -i ${id} -1 ${read1} -2 ${read2}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(prepare_fastqs.py -v)
        container: ${task.container}
    END_VERSIONS
    """

}
