process ANNOTATE_SV {
    tag "$meta.id"
    label 'process_low'
    container "ghcr.io/dhslab/docker-vep:release_105"

    input:
    tuple val(meta), path(align_output)
    tuple val(chromoseq_inputs), path("*")
    val(chromoseq_parameters)

    output:
    tuple val(meta), path("${meta.id}.sv_vep.vcf.gz*", arity:'2'), emit: vcf
    path "versions.yml",    emit: versions

    script:
    """
    /usr/bin/perl -I /opt/lib/perl/VEP/Plugins /opt/vep/src/ensembl-vep/vep --format vcf --vcf \\
    --fasta ${params.genome} --flag_pick --symbol --distance ${chromoseq_parameters.sv_annotation_distance} \\
    --term SO -o ${meta.id}.sv_vep.vcf -i ${meta.id}.sv.vcf.gz \\
    --custom ${chromoseq_inputs.cytobands},cytobands,bed --custom ${chromoseq_inputs.sv_bed},KnownSvGenes,bed \\
    --offline --cache --dir ${chromoseq_inputs.vepcache} && \\
    bgzip -c ${meta.id}.sv_vep.vcf > ${meta.id}.sv_vep.vcf.gz && \\
    tabix -p vcf ${meta.id}.sv_vep.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep --version)
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}.sv_vep.vcf.gz"
    touch "${meta.id}.sv_vep.vcf.gz.tbi"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vep: \$(/opt/vep/src/ensembl-vep/vep --version)
    END_VERSIONS
    """
}
