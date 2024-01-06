process COMBINE_RUNINFO {
    label 'process_tiny'
    container "ghcr.io/dhslab/docker-python3:231224"
    
    input:
    tuple val(id), path('runinfo?.csv')

    output:
    tuple val(id), path('runinfo.json'), emit: runinfo
    path 'versions.yml', emit: versions

    script: 
    """
    csv2json.py -i ${id} runinfo*.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csv2json.py: \$(csv2json.py -v)
    END_VERSIONS
    """

    stub:
    """
    csv2json.py -i ${id} -u runinfo*.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csv2json.py: \$(csv2json.py -v)
    END_VERSIONS
    """

}
