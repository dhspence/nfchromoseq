//
// This file holds several functions specific to the workflow/nfchromoseq.nf in the dhslab/nfchromoseq pipeline
//

@Grab(group='org.apache.commons', module='commons-csv', version='1.8')
import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine
import org.apache.commons.csv.CSVFormat
import org.apache.commons.csv.CSVParser
import org.apache.commons.csv.CSVRecord
import java.nio.file.Files
import java.nio.file.Paths

class WorkflowNfchromoseq {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        genomeExistsError(params, log)
        if (!params.fasta) {
            log.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
            System.exit(1)
        }

        // Check that the mastersheet exists
        if (!params.outdir || !Files.exists(Paths.get(params.outdir))) {
            log.error "Output directory not specified."
            System.exit(1)
        }

        // Check that the mastersheet exists
        if (!params.input || params.input==true || !Files.exists(Paths.get(params.input))) {
            log.error "Mastersheet: (${params.input}) not specified with e.g. '--input mastersheet.csv' or via a detectable config file."
            System.exit(1)
        }

        // run validateMastersheet
        if (validateMastersheet(params,log) == false) {
            log.error "Mastersheet validation failed."
            System.exit(1)
        }

        // validate all files in params.dragen_inputs
        params.dragen_inputs.each { key, value ->
            if (value == null || !Files.exists(Paths.get(value))) {
                log.error "Invalid path in 'dragen_inputs' field: ${key}"
                System.exit(1)
            }
        }

        // validate all files in params.chromoseq_inputs
        params.chromoseq_inputs.each { key, value ->
            if (value == null || !Files.exists(Paths.get(value))) {
                log.error "Invalid path in 'chromoseq_inputs' field: ${key}"
                System.exit(1)
            }
        }

    }

    //
    // Validate mastersheet
    //
    public static boolean validateMastersheet(params,log) {
        def csvFile = new File(params.input)
        // Check if the file exists
        if (!csvFile.exists()) {
            log.error "Workflow validation error: Mastersheet not found: ${csvFile.path}"
            System.exit(1)
        }

        // Read the file and check headers
        csvFile.withReader { reader ->
            CSVParser parser = CSVFormat.DEFAULT.withFirstRecordAsHeader().parse(reader)
            def headers = parser.headerMap.keySet()

            if (!headers.contains('id')) {
                log.error "Workflow validation error: No 'id' field found in mastersheet."
                System.exit(1)
            }

            for (CSVRecord record : parser) {

                def id = record.get('id')
                if (id.trim().isEmpty()) {
                    log.error "Workflow validation error: Invalid 'id' field: cannot be null or empty."
                    System.exit(1)
                }
                // Validate 'id' field
                if (id.contains(" ")) {
                    log.error "Mastersheet error: id field should not contain spaces."
                    System.exit(1)
                }
                
                // Validate 'lanes' field
                if (headers.contains('lanes') || headers.contains('index')){
                    
                    if ((record.get('lanes').trim() && !record.get('index').trim()) || (!record.get('lanes').trim() && record.get('index').trim())){
                        log.error "Workflow validation error: index and lanes must both be specified."
                        System.exit(1)
                    }

                    if (record.get('lanes').trim() && record.get('index').trim() && params.rundir == null){
                        log.error "Workflow validation error: rundir must be given if lanes and index are specified."
                        System.exit(1)
                    }

                    if (!record.get('lanes').matches("\\d+(,\\d+)*")) {
                        log.error "Workflow validation error: Invalid format in 'lanes' mastersheet field (e.g., 1,2,3)."
                        System.exit(1)
                    }

                    // Validate 'index' field
                    if (headers.contains('index')){
                        def index = record.get('index')
                        if (index.find("[^ACGT\\-]") || index.find(" ")) {
                            log.error "Workflow validation error: Invalid format in 'index': ${index}."
                            System.exit(1)
                        }
                    }
                }

                headers.eachWithIndex { header, index ->
                    if (header in ['fastq_list','demux_path', 'dragen_path', 'read1', 'read2']){
                        def filepath = record.get(header).trim()
                        if (!Files.exists(Paths.get(filepath))) {
                            log.error "Workflow validation error: Invalid path in '$header' field: ${filepath}"
                            System.exit(1)
                        }
                    }
                }

                // Check at least one of demux_path, dragen_dir, or index are valid
                if (!['demux_path', 'dragen_path', 'fastq_list', 'index','read1','read2'].any { headers.contains(it) && record.get(it).trim() }) {
                    log.error "Workflow validation error: At least one of 'demux_path', 'dragen_path', 'index', or 'read1' and 'read2' must be valid."
                    System.exit(1)
                }

            }
        }
        return params
    }

    
    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += "    </dl>\n"
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/','-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += "data: |\n"
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Generate methods description for MultiQC
    //

    public static String toolCitationText(params) {

        // TODO nf-core: Optionally add in-text citation tools to this list.
        // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
        // Uncomment function in methodsDescriptionText to render in MultiQC report
        def citation_text = [
                "Tools used in the workflow included:",
                "FastQC (Andrews 2010),",
                "MultiQC (Ewels et al. 2016)",
                "."
            ].join(' ').trim()

        return citation_text
    }

    public static String toolBibliographyText(params) {

        // TODO Optionally add bibliographic entries to this list.
        // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
        // Uncomment function in methodsDescriptionText to render in MultiQC report
        def reference_text = [
                "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
                "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
            ].join(' ').trim()

        return reference_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml, params) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta["manifest_map"] = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

        // Tool references
        meta["tool_citations"] = ""
        meta["tool_bibliography"] = ""

        // TODO Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
        //meta["tool_citations"] = toolCitationText(params).replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
        //meta["tool_bibliography"] = toolBibliographyText(params)


        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }
}
