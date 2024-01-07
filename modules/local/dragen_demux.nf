process DRAGEN_DEMUX {
    label 'dragen'
    container "${params.dragen_container}"

    input:
    path samplesheet
    path rundir
    path demuxdir

    output:
    path ('fastq_list.csv'), emit: fastqlist
    path "versions.yml",    emit: versions

    script:
    def first_tile = params.bcl_first_tile ? " --first-tile-only true" : ""
    """
    /opt/edico/bin/dragen --bcl-conversion-only true --bcl-only-matched-reads true --strict-mode true${first_tile} \\
    --sample-sheet ${samplesheet} --bcl-input-directory ${rundir} \\
    --output-directory \$(realpath ${demuxdir}) > ./demux_log.txt && \\
    cp demux_log.txt ${demuxdir}/Reports/ && \\
    cp ${rundir}/RunParameters.xml ${demuxdir}/Reports/ && \\
    cp ${demuxdir}/Reports/fastq_list.csv fastq_list.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(/opt/edico/bin/dragen --version | tail -n 1 | cut -d ' ' -f 3)
    END_VERSIONS
    """

    stub:
    def first_tile = params.bcl_first_tile ? " --first-tile-only true" : ""
    """
    echo ${first_tile} && \\
    mkdir \$(realpath ${demuxdir}) && \\
    mkdir ${demuxdir}/Reports/ && \\
    cp ${projectDir}/assets/stub/demux_fastq/Reports/fastq_list.csv ${demuxdir}/Reports/fastq_list.csv && \
    cp ${demuxdir}/Reports/fastq_list.csv .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen: \$(cat ${projectDir}/assets/stub/versions/dragen_version.txt)
    END_VERSIONS
    """

}