process ANNOTATE_VARIANTS {
    tag "$meta.id"
    label 'process_low'
    label 'final_output'
    container "ghcr.io/dhslab/docker-vep:release_105"

    input:
    tuple val(meta), path(vcf)
    tuple val(chromoseq_inputs), path("*")
    val(chromoseq_parameters)

    output:
    tuple val(meta), path("${meta.id}.annotated.vcf.gz*", arity: '2'), emit: vcf
    path "versions.yml",    emit: versions

    script:
    """
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /opt/vep/src/ensembl-vep/vep \\
    --format vcf --vcf --fasta ${params.fasta} --hgvs --symbol --term SO --flag_pick -o ${meta.id}.annotated.vcf \\
    -i ${meta.id}.smallvariants.vcf.gz --custom ${chromoseq_inputs.cytobands},cytobands,bed --custom ${chromoseq_inputs.custom_annotation_vcf},${chromoseq_parameters.custom_annotation_parameters} --offline --cache --max_af --dir ${chromoseq_inputs.vepcache} && \\
    bgzip -c ${meta.id}.annotated.vcf > ${meta.id}.annotated.vcf.gz && \\
    tabix -p vcf ${meta.id}.annotated.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep 2>&1 | grep ensembl-vep | cut -d ':' -f 2 | sed 's/\s*//g')
    END_VERSIONS
    /opt/vep/src/ensembl-vep/vep --dir ${chromoseq_inputs.vepcache} --show_cache_info | awk '{ print "    "\$1": "\$2; }' >> versions.yml
    """

    stub:
    """
    touch "${meta.id}.annotated.vcf.gz"
    touch "${meta.id}.annotated.vcf.gz.tbi"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
    \$(cat $projectDir/assets/stub/versions/vep_version.yaml)
    END_VERSIONS
    """
}
