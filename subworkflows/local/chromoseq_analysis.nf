include { GET_COVERAGE          } from '../../modules/local/get_coverage.nf'
include { RUN_ICHOR             } from '../../modules/local/run_ichor.nf'
include { ANNOTATE_VARIANTS     } from '../../modules/local/annotate_variants.nf'
include { FILTER_VARIANTS       } from '../../modules/local/filter_variants.nf'
include { FILTER_CNV            } from '../../modules/local/filter_cnv.nf'
include { ANNOTATE_SV           } from '../../modules/local/annotate_sv.nf'
include { FILTER_SV             } from '../../modules/local/filter_sv.nf'
include { RUN_HAPLOTECT         } from '../../modules/local/run_haplotect.nf'
include { MAKE_REPORT           } from '../../modules/local/make_report.nf'

// function to get the sample sex from the mapping_metrics.csv or cnv_metrics.csv file
def getSampleSex(List<String> filePaths) {
    result = null
    // Iterate through the list of file paths
    filePaths.each { filePath ->
        def myfile = file(filePath)
        // Check if the file is 'mapping_metrics.csv'
        if (myfile.name.endsWith('.mapping_metrics.csv') && myfile.exists()) {
            // Read the file line by line
            myfile.eachLine { line ->
                def columns = line.split(',') // Assuming comma-separated values
                if (columns[2] == 'Provided sex chromosome ploidy' && (columns[3] == 'XX' || columns[3] == 'XY')) {
                    result = (columns[3] == 'XX') ? 'female' : 'male'
                }
            }
        }
    }

    filePaths.each { filePath ->
        def myfile = file(filePath)

        // Check if the file is 'mapping_metrics.csv'
        if (myfile.name.endsWith('.cnv_metrics.csv') && myfile.exists()) {
            // Read the file line by line
            myfile.eachLine { line ->
                def columns = line.split(',') // Assuming comma-separated values
                if (columns[0] == 'SEX GENOTYPER' && (columns[3] == 'XX' || columns[3] == 'XY')) {
                    result = (columns[3] == 'XX') ? 'female' : 'male'
                }
            }
        }
    }
    return result
}

workflow CHROMOSEQ_ANALYSIS {
    take:
    dragen_outputs // channel.of(<dragen output paths>)
    chromoseq_inputs // channel.of(<staged chromoseq input files>)   
    chromoseq_parameters // <chromoseq parameters>

    main:
    ch_versions = Channel.empty()
    
    // get sex from the dragen output if its not provided in the mastersheet
    dragen_outputs
    .map { meta, outputs -> 
        if (meta.sex == null){
            meta.sex = getSampleSex(outputs)
        }
        [ meta, outputs]
    }
    .set { ch_dragen_outputs }

    GET_COVERAGE ( 
        ch_dragen_outputs,
        chromoseq_inputs,
        chromoseq_parameters
    )
    ch_versions.mix(GET_COVERAGE.out.versions)

    RUN_ICHOR ( 
        ch_dragen_outputs,
        chromoseq_inputs,
        chromoseq_parameters
    )
 
    FILTER_VARIANTS (
        ch_dragen_outputs,
        chromoseq_inputs,
        chromoseq_parameters
    )
    ch_versions.mix(FILTER_VARIANTS.out.versions)

    ch_filtered_variants_vcf = FILTER_VARIANTS.out.vcf
    ANNOTATE_VARIANTS (
        ch_filtered_variants_vcf,
        chromoseq_inputs,
        chromoseq_parameters
    )
    ch_versions.mix(ANNOTATE_VARIANTS.out.versions)

    ch_cnv_vcf = RUN_ICHOR.out.seg.join(ch_dragen_outputs)
    FILTER_CNV (
        ch_cnv_vcf,
        chromoseq_inputs,
        chromoseq_parameters
    )
    ch_versions.mix(FILTER_CNV.out.versions)

    ANNOTATE_SV (
        ch_dragen_outputs,
        chromoseq_inputs,
        chromoseq_parameters
    )
    ch_versions.mix(ANNOTATE_SV.out.versions)

    ch_annotated_sv_vcf = ANNOTATE_SV.out.vcf
    FILTER_SV (
        ch_annotated_sv_vcf,
        chromoseq_inputs,
        chromoseq_parameters
    )
    ch_versions.mix(FILTER_SV.out.versions)

    RUN_HAPLOTECT (
        ch_dragen_outputs,
        chromoseq_inputs
    )

    make_report_input = RUN_ICHOR.out.seg
    .concat( GET_COVERAGE.out.report, ANNOTATE_VARIANTS.out.vcf, FILTER_CNV.out.vcf, FILTER_SV.out.vcf, RUN_HAPLOTECT.out.out_file, ch_dragen_outputs)
    .transpose()
    .groupTuple()

    MAKE_REPORT (
        make_report_input, 
        chromoseq_inputs, 
        chromoseq_parameters
    )
    ch_versions.mix(MAKE_REPORT.out.versions)

    emit:
    report = MAKE_REPORT.out.json_report
    versions = ch_versions

}