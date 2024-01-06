workflow INPUT_CHECK {
    take:
        samplesheet // file: /path/to/samplesheet.csv

    main:
        SAMPLESHEET_CHECK ( samplesheet )
            .csv
            .splitCsv ( header:true, sep:',', quote:'"' )
            .map { create_samplesheet(it) }
            .set { meta }

    emit:
        meta                                     // channel: [ val(meta), [ reads ] ]
        versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}

process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    label 'process_tiny'

    input:
    path samplesheet

    output:
    path('samplesheet.valid.csv'), emit: csv
    path('versions.yml'),          emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    check_samplesheet.py ${samplesheet} samplesheet.valid.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
    
}

def create_samplesheet(LinkedHashMap row) {
    // create meta map
    def meta = [:]
    meta.id             = row.id ?: null
    meta.lanes          = row.lanes ? row.lanes.split(',').collect{ it as int } : null

    def i7 = null
    def i5 = null
    if (row.index && row.index.contains('-')) {
        (i7, i5) = row.index.split('-', 2) // Splitting with limit to handle unexpected formats
    }
    meta.i7index        = i7
    meta.i5index        = i5
    meta.exceptions     = row.exceptions ?: null
    meta.mrn            = row.mrn ?: null
    meta.accession      = row.accession ?: null  
    meta.specimen       = row.specimen ?: null
    meta.dob            = row.dob ?: null
    meta.sex            = row.sex ?: null
    meta.fastq_list     = row.fastq_list ?: null
    meta.demux_path     = row.demux_path ?: null
    meta.dragen_path    = row.dragen_path ?: null
    meta.read1          = row.read1 ?: null
    meta.read2          = row.read2 ?: null

    return meta
}
