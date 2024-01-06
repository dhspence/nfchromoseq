workflow DRAGEN_ALIGN {
    take:
    ch_fastq_inputs
    ch_cram_inputs
    ch_dragen_inputs
    chromoseq_parameters
        
    main:
    ch_versions = Channel.empty()
    ch_dragen_paths = Channel.empty()

        // make demux samplesheet
    DRAGEN_ALIGN_CRAM(ch_cram_inputs, ch_dragen_inputs, chromoseq_parameters)
    ch_dragen_paths = ch_dragen_paths.mix(DRAGEN_ALIGN_CRAM.out.dragen_output)        
    ch_versions = ch_versions.mix(DRAGEN_ALIGN_CRAM.out.versions)

    DRAGEN_ALIGN_FASTQLIST(ch_fastq_inputs, ch_dragen_inputs, chromoseq_parameters)
    ch_dragen_paths = ch_dragen_paths.mix(DRAGEN_ALIGN_FASTQLIST.out.dragen_output)        
    ch_versions = ch_versions.mix(DRAGEN_ALIGN_FASTQLIST.out.versions)

    emit:
    dragen_output = ch_dragen_paths
    versions = ch_versions
}

process DRAGEN_ALIGN_CRAM {
    label 'dragen'
    container "$params.dragen_container"
    publishDir "$params.outdir/${meta.id}/dragen", mode:'copy'

    input:
        tuple val(meta), path(cram)
        tuple val(dragen_inputs), path("*")
        val(chromoseq_parameters)

    output:
        tuple val(meta), path ("${meta.id}*"), emit: dragen_output
        path "versions.yml"                  , emit: versions

    script:
        def specified_sex = meta.sex != null ? "--sample-sex ${meta.sex}" : ""
        """
        /opt/edico/bin/dragen -r ${dragen_inputs.reference} ${specified_sex} \\
        --read-trimmers adapter --trim-adapter-read1 ${dragen_inputs.dragen_adapter1} --trim-adapter-read2 ${dragen_inputs.dragen_adapter2} \\
        --tumor-cram-input ${cram} --cram-reference ${params.input_cram_reference} \\
        --qc-coverage-region-1 ${dragen_inputs.cna_windows} --qc-coverage-filters-1 'mapq<0,bq<0' \\
        --qc-coverage-ignore-overlaps=true \\
        --gc-metrics-enable true --enable-map-align-output true --enable-bam-indexing true --enable-duplicate-marking true \\
        --enable-variant-caller true --dbsnp ${dragen_inputs.dbsnp} --vc-somatic-hotspots ${dragen_inputs.hotspots} --vc-systematic-noise ${dragen_inputs.snv_noisefile} --vc-enable-triallelic-filter false --vc-combine-phased-variants-distance 3 \\
        --enable-cnv true --cnv-somatic-enable-het-calling true --cnv-enable-ref-calls false --cnv-merge-distance ${chromoseq_parameters.dragen_cnv_merge_distance} --cnv-filter-length ${chromoseq_parameters.dragen_cnv_filter_length} --cnv-population-b-allele-vcf ${dragen_inputs.pop_af_vcf} \\
        --enable-sv true --sv-output-contigs true --sv-hyper-sensitivity true --sv-min-edge-observations 2 --sv-min-candidate-spanning-count 1 \\
        --sv-use-overlap-pair-evidence true --sv-enable-somatic-ins-tandup-hotspot-regions true --sv-systematic-noise ${dragen_inputs.sv_noisefile} \\
        --output-format CRAM --output-directory ./ --force --output-file-prefix ${meta.id} --intermediate-results-dir ${params.dragen_staging_path} &> ./${meta.id}_dragen.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            dragen: \$(/opt/edico/bin/dragen --version | tail -n 1 | cut -d ' ' -f 3)
        END_VERSIONS
        """

    stub:
        """
        cp ${projectDir}/assets/stub/${meta.id}/* .
        echo -e "${task.process}:\n\tdragen: "\$(cat $projectDir/assets/stub/versions/dragen_version.txt) > versions.yml
        """

}

process DRAGEN_ALIGN_FASTQLIST {
    label 'dragen'
    container "$params.dragen_container"
    publishDir "$params.outdir/${meta.id}/dragen", mode:'copy'

    input:
        tuple val(meta), val(fastqlist), path(read1), path(read2), path(runinfo)
        tuple val(dragen_inputs), path("*")
        val(chromoseq_parameters)

    output:
        tuple val(meta), path ("${meta.id}*"), emit: dragen_output
        path "versions.yml"                  , emit: versions

    script:
        def specified_sex = meta.sex != null ? "--sample-sex ${meta.sex}" : ""
        """
        echo "${fastqlist}" > fastqlist.csv && \\
        cp ${runinfo} ${meta.id}.run_info.json && \\
        /opt/edico/bin/dragen -r ${dragen_inputs.reference} ${specified_sex}\\
        --read-trimmers adapter --trim-adapter-read1 ${dragen_inputs.dragen_adapter1} --trim-adapter-read2 ${dragen_inputs.dragen_adapter2} \\
        --tumor-fastq-list fastqlist.csv --tumor-fastq-list-sample-id ${meta.id} \\
        --qc-coverage-region-1 ${dragen_inputs.cna_windows} --qc-coverage-filters-1 'mapq<0,bq<0' \\
        --qc-coverage-ignore-overlaps=true \\
        --gc-metrics-enable true --enable-map-align-output true --enable-bam-indexing true --enable-duplicate-marking true \\
        --enable-variant-caller true --dbsnp ${dragen_inputs.dbsnp} --vc-somatic-hotspots ${dragen_inputs.hotspots} --vc-systematic-noise ${dragen_inputs.snv_noisefile} --vc-enable-triallelic-filter false --vc-combine-phased-variants-distance 3 \\
        --enable-cnv true --cnv-somatic-enable-het-calling true --cnv-enable-ref-calls false --cnv-merge-distance ${chromoseq_parameters.dragen_cnv_merge_distance} --cnv-filter-length ${chromoseq_parameters.dragen_cnv_filter_length} --cnv-population-b-allele-vcf ${dragen_inputs.pop_af_vcf} \\
        --enable-sv true --sv-output-contigs true --sv-hyper-sensitivity true --sv-min-edge-observations 2 --sv-min-candidate-spanning-count 1 \\
        --sv-use-overlap-pair-evidence true --sv-enable-somatic-ins-tandup-hotspot-regions true --sv-systematic-noise ${dragen_inputs.sv_noisefile} \\
        --output-format CRAM --output-directory ./ --force --output-file-prefix ${meta.id} --intermediate-results-dir ${params.dragen_staging_path} &> ./${meta.id}_dragen.log

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            dragen: \$(/opt/edico/bin/dragen --version | tail -n 1 | cut -d ' ' -f 3)
        END_VERSIONS
        """

    stub:
        """
        cp $projectDir/assets/stub/dragen_path/${meta.id}* .
        echo -e "${task.process}:\n\tdragen: "\$(cat $projectDir/assets/stub/versions/dragen_version.txt) > versions.yml
        """

}
