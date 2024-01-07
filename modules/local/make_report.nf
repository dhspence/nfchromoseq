process MAKE_REPORT { 
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-cleutils:2"
    
    publishDir = [
        path: { "${params.outdir}/${meta.id}" }, 
        mode: params.publish_dir_mode, 
        pattern: "*.chromoseq.{txt,json}"
    ]

    input:
    tuple val(meta), path(makereportinput)
    tuple val(chromoseq_inputs), path("*")
    val(chromoseq_parameters)

    output:
    path("${meta.id}.chromoseq.txt")            , emit: text_report
    path("${meta.id}.chromoseq.json")           , emit: json_report
    path("versions.yml")                        , emit: versions

    script:
    """
    make_cs_report.py -p ${chromoseq_parameters.max_pop_af} --sex ${meta.sex} --DOB ${meta.dob} --exception "${meta.exceptions}" \\
    --genebed ${chromoseq_inputs.gene_bed} --svbed ${chromoseq_inputs.sv_bed} --svtargets ${chromoseq_inputs.sv_targets} --qc ${chromoseq_inputs.assay_specifications} -n ${meta.id} -d ./

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(make_cs_report.py --version)
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}.chromoseq.txt"
    touch "${meta.id}.chromoseq.json"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(make_cs_report.py --version)
    END_VERSIONS
    """
}