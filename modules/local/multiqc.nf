process MULTIQC {
    label 'process_low'
    container "quay.io/biocontainers/multiqc:1.13--pyhdfd78af_0"
    // publishDir "$params.outdir", mode:'copy'

    input:
    path  multiqc_files, stageAs: "?/*"

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data

    script:
    """
    multiqc .
    """

}