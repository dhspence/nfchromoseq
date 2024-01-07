/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowNfchromoseq.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK            } from '../subworkflows/local/input_check.nf'
include { CHROMOSEQ_ANALYSIS     } from '../subworkflows/local/chromoseq_analysis.nf'
include { GATHER_FASTQS          } from '../subworkflows/local/gather_fastqs.nf'
include { DEMUX                  } from '../subworkflows/local/demux.nf'
include { DRAGEN_ALIGN           } from '../subworkflows/local/dragen_align.nf'
include { COMBINE_RUNINFO        } from '../modules/local/combine_runinfo.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CUSTOM HELPER FUNCTIONS FOR MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// This function 'stages' a set of files defined by a map of key:filepath pairs.
// It returns a tuple: a map of key:filename pairs and list of file paths.
// This can be used to generate a value Channel that can be used as input to a process
// that accepts a tuple val(map), path("*") so map.key refers to the appropriate linked file.
def stageFileset(Map filePathMap) {
    def basePathMap = [:]
    def filePathsList = []

    filePathMap.each { key, value ->

        def filepath = file(value, checkExists: true)
        if (filepath.exists()) {
            // Add basename and key to the map
            basePathMap[key] = value.split('/')[-1]
            // Add file path to the list
            filePathsList << filepath
        } else {
            println "Warning: File at '${value}' for key '${key}' does not exist."
        }
    }
    return [basePathMap, filePathsList]
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow NFCHROMOSEQ {

    ch_versions = Channel.empty()

    // Mastersheet channel
    ch_mastersheet = Channel.fromPath(params.input)

    // stage chromoseq-specific inputs
    ch_chromoseq_inputs = Channel.value(stageFileset(params.chromoseq_inputs))

    // add cna windows and hotspot vcf to params.dragen_inputs and stage dragen input files
    params.dragen_inputs.hotspots = params.chromoseq_inputs.hotspots
    params.dragen_inputs.hotspots_index = params.chromoseq_inputs.hotspots_index
    params.dragen_inputs.cna_windows = params.chromoseq_inputs.cna_windows
    ch_dragen_inputs = Channel.value(stageFileset(params.dragen_inputs))

    // holds fastq list files for processing and alignment
    ch_fastqs = Channel.empty()

    // holds fastq list files for processing and alignment
    ch_cramstorealign = Channel.empty()

    // holds fastq lists and ids for alignment
    ch_align_inputs = Channel.empty()

    // holds dragen output paths
    ch_dragen_output = Channel.empty()

    // holds processed dragen output paths and metadata for analysiss
    ch_analysis_inputs = Channel.empty()    

    // check the samplesheet. This just gets the meta values.
    INPUT_CHECK(ch_mastersheet)
    ch_sample_meta = INPUT_CHECK.out.meta
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // only run demux if a rundir is passed
    if (params.rundir && params.run_demux == true){
        // do demux
        DEMUX(ch_mastersheet, Channel.fromPath(params.rundir), Channel.fromPath(params.demuxdir))
        ch_fastqs = ch_fastqs.mix(DEMUX.out.fastqlist)
        ch_versions = ch_versions.mix(DEMUX.out.versions)
    }    

    // get fastq lists from mastersheet, including from fastq_list.csv, demux_directories, and crams.
    // FASTQs are generated from all of these sources and compiled into a single list for dragen align.
    // Runinfo data are also compiled.
    GATHER_FASTQS(ch_sample_meta)
    ch_fastqs = ch_fastqs.mix(GATHER_FASTQS.out.fastqlist)
    ch_versions = ch_versions.mix(GATHER_FASTQS.out.versions)

    // now combine all fastqs into a single object with the runinfo
    ch_fastqs
    .groupTuple()
    .map { id, fqlist, runinfo -> 
            def fileList = ['RGID,RGSM,RGLB,Lane,Read1File,Read2File']
            def read1 = []
            def read2 = []

            // Create data rows
            for (int i = 0; i < fqlist.size(); i++) {
                def row = fqlist[i]
                read1 << file(row[3])
                read2 << file(row[4])
                fileList << [ row[0], id, row[1], row[2], row[3].toString().split('/')[-1], row[4].toString().split('/')[-1] ].join(',')
            }
            return [ id, fileList.join('\n'), read1, read2, runinfo ]
    }
    .set { ch_fastqs }

    // the run info needs to be consolidate, so first we have to separate it out then recombine it
    ch_fastqs
    .map { id, list, read1, read2, runinfo -> [ id, runinfo ] }
    .set { ch_runinfo }

    COMBINE_RUNINFO(ch_runinfo)
    ch_versions = ch_versions.mix(COMBINE_RUNINFO.out.versions)

    ch_fastqs
    .map { id, list, read1, read2, runinfo -> [ id, list, read1, read2 ] }
    .join(COMBINE_RUNINFO.out.runinfo)
    .set { ch_fastqs }

    // now add the full metadata
    ch_sample_meta
    .map { meta -> 
        new_meta = meta.subMap('id', 'mrn', 'accession', 'specimen', 'dob', 'sex','exceptions')
        [meta.id, new_meta]
    }
    .unique()
    .join(ch_fastqs)
    .map {id, meta, list, read1, read2, runinfo -> [meta, list, read1, read2, runinfo]}
    .set { ch_align_input }

    // finally, get crams that need to be realigned, if realignment was specified
    ch_sample_meta
    .map { meta -> 
        if (params.realign && !meta.i5index && !meta.i7index && !meta.fastq_list && !meta.read1 && !meta.read2) {
            new_meta = meta.subMap('id', 'mrn', 'accession', 'specimen', 'dob', 'sex','exceptions')
            if (meta.dragen_path) {
                if (file(meta.dragen_path).isDirectory()) {
                    def cram = file("${meta.dragen_path}/${meta.id}_tumor.cram", checkExists: true)
                    [ new_meta, cram ]
                } else if (file(meta.dragen_path).isFile()) {
                    [ new_meta, file(meta.dragen_path) ]
                }
            }
        }
    }
    .set { ch_cramtorealign }

    if (params.debug){
        ch_align_input.view()
    }

    if (params.run_align == true) {
        DRAGEN_ALIGN(ch_align_input, ch_cramtorealign, ch_dragen_inputs, params.chromoseq_parameters)
        ch_dragen_output = DRAGEN_ALIGN.out.dragen_output
        ch_versions = ch_versions.mix(DRAGEN_ALIGN.out.versions)
    }

    // add any dragen paths to analysis list, if realign was not specified
    if (!params.realign){
        ch_sample_meta
        .map { meta -> 
            if (!meta.read1 && !meta.read2 && !meta.fastq_list && meta.dragen_path != null){
                def files = file("${meta.dragen_path}/${meta.id}*")
                return([ meta, files ])
            }
        }
        .concat(ch_dragen_output)
        .set { ch_dragen_output }
    }

    if (params.debug){
        ch_dragen_output.view()
    }

    // run analysis
    if (params.run_analysis == true) {
        CHROMOSEQ_ANALYSIS(ch_dragen_output, ch_chromoseq_inputs, params.chromoseq_parameters)
        ch_versions = ch_versions.mix(CHROMOSEQ_ANALYSIS.out.versions)
    }

//    ch_qc_reports = CHROMOSEQ_ANALYSIS.out.mqc_report

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowNfchromoseq.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowNfchromoseq.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
//    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
