process FASTQ_FROM_CRAM {
    label 'process_high_multithread'
    container "ghcr.io/dhslab/docker-htslib:231224"

    input:
        tuple val(id), path(cramfile)

    output:
        tuple val(id), path("*.R1.fastq.gz"), path("*.R2.fastq.gz"), emit: fastqs
        path("versions.yml"),    emit: versions

    script:
        def threads = (task.cpus / 2) - 1

        """
        samtools collate -@ ${threads} --reference ${params.input_cram_reference} -u -O ${cramfile} | \
        samtools fastq -@ ${threads} -1 \$(basename ${cramfile} .cram).R1.fastq.gz -2 \$(basename ${cramfile} .cram).R2.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)
        END_VERSIONS
        """

    stub:
        """
        samtools collate --reference ${params.input_cram_reference} -u -O ${projectDir}/assets/stub/${id}/*.cram | \
        samtools fastq -1 \$(basename ${cramfile} .cram).R1.fastq.gz -2 \$(basename ${cramfile} .cram).R2.fastq.gz

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            samtools: \$(samtools --version | head -n 1 | cut -d ' ' -f 2)
        END_VERSIONS
    """

}
