process RUN_ICHOR {
    tag "$meta.id"
    label 'process_medium'
    label 'final_output'
    container "mgibio/chromoseq:v1.1"

    input:
    tuple val(meta), path(align_output)
    tuple val(chromoseq_inputs), path("*")
    val(chromoseq_parameters)
    
    output:
    tuple val(meta), path ("*.segs.txt")                  , emit: seg
    tuple val(meta), path ("*_genomeWide_all_sols.pdf")   , emit: allgenomewide_pdf
    tuple val(meta), path ("*_genomeWide.pdf")            , emit: genomewide_pdf
    tuple val(meta), path ("*_genomeWideCorrection.pdf")  , emit: genomewideCorrection_pdf
 
    script:
    """
    set -eo pipefail && \\
    /bin/cat ${meta.id}.qc-coverage-region-1_read_cov_report.bed | awk 'NR>1' | sort -k 1V,1 -k 2n,2 | \\
    awk -v window=${chromoseq_parameters.ichor_bin_size} 'BEGIN { chr=""; } { if (\$1!=chr){ printf("fixedStep chrom=%s start=1 step=%d span=%d\\n",\$1,window,window); chr=\$1; } print \$4; }' > "${meta.id}.tumor.wig" && \\
    /usr/local/bin/Rscript /usr/local/bin/ichorCNA/scripts/runIchorCNA.R --id ${meta.id} \\
    --WIG "${meta.id}.tumor.wig" --ploidy "c(${chromoseq_parameters.ichor_ploidy_levels})" --normal "c(0.1,0.5,.85)" --maxCN 3 \\
    --gcWig ${chromoseq_inputs.gc} \\
    --mapWig ${chromoseq_inputs.mappability} \\
    --centromere ${chromoseq_inputs.centromeres} \\
    --normalPanel ${chromoseq_inputs.pon} \\
    --genomeBuild hg38 \\
    --sex ${meta.sex} \\
    --includeHOMD False --chrs 'c(1:22, \"X\", \"Y\")' --chrTrain 'c(1:22)' --fracReadsInChrYForMale 0.0005 \\
    --estimateNormal True --estimatePloidy True --estimateScPrevalence True \\
    --txnE 0.999999 --txnStrength 1000000 --genomeStyle UCSC --outDir ./ --libdir /usr/local/bin/ichorCNA/ && \\
    awk -v G=${meta.sex} '\$2!~/Y/ || G=="male"' "${meta.id}.seg.txt" > "${meta.id}.segs.txt" && \\
    mv ${meta.id}/*.pdf .
    """

    stub:
    """
    touch "${meta.id}.segs.txt"
    touch "${meta.id}_genomeWide_all_sols.pdf"
    touch "${meta.id}_genomeWideCorrection.pdf"
    touch "${meta.id}_genomeWide.pdf"
    """

}