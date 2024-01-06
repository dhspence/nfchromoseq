process SPLIT_CRAM {
    label 'process_medium_multithread'
    container "ghcr.io/dhslab/docker-htslib:231224"

    input:
        tuple val(id), path(cramfile)
    
    output:
        tuple val(id), path("*.cram"),      emit: crams
        path("versions.yml"),               emit: versions

    script:
        """
        samtools split -@ ${task.cpus} -f '%!.cram' --write-index --reference ${params.input_cram_reference} ${cramfile}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)
        """

    stub:
        """
        cp ${projectDir}/assets/stub/fastq_from_cram/split_crams/*.cram .
        cp ${projectDir}/assets/stub/fastq_from_cram/split_crams/*.crai .

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)
        END_VERSIONS
        """
}
