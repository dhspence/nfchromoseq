process RUN_HAPLOTECT {
    tag "$meta.id"
    label 'process_medium'
    label 'final_output'
    container "ghcr.io/dhslab/docker-haplotect:231129"

    input:
    tuple val(meta), path(align_output)
    tuple val(chromoseq_inputs), path("*")

    output:
    tuple val(meta), path("${meta.id}.haplotect.txt"), emit: out_file

    script:
    """
    haplotect -e ${meta.id}"_tumor.cram" ${chromoseq_inputs.haplotect} ${params.fasta} > ${meta.id}".haplotect.txt"
    """

    stub:
    """
    touch "${meta.id}.haplotect.txt"
    """
}
