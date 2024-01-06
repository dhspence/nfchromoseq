include { MAKE_FASTQLIST      } from '../../modules/local/make_fastqlist.nf'
include { GET_FASTQ_INFO      } from '../../modules/local/get_fastq_info.nf'
include { SPLIT_CRAM          } from '../../modules/local/split_cram.nf'
include { FASTQ_FROM_CRAM     } from '../../modules/local/fastq_from_cram.nf'

workflow GATHER_FASTQS {
    take:
    ch_sample_meta

    main:

    ch_fastqs = Channel.empty()
    ch_runinfo = Channel.empty()
    ch_versions = Channel.empty()

        // get fastq_lists from demux dir
    ch_sample_meta
    .map { meta -> 
        if (meta.demux_path != null && file(meta.demux_path).isDirectory()) {
            [ meta.id, meta.demux_path ] 
        }
    }
    .set { ch_demuxpaths }

    // combine with fastq_lists from samplesheet
    ch_sample_meta
    .map { meta -> 
        if (meta.fastq_list && file(meta.fastq_list).isFile()) {
            [ meta.id, meta.fastq_list ] 
        }
    }
    .concat(ch_demuxpaths)
    .set { ch_fastqlists }    

    // get fastqs from read1, read2 in mastersheet
    ch_sample_meta
    .map { meta -> 
        if (meta.read1 && file(meta.read1).isFile() && meta.read2 && file(meta.read2).isFile()) {
            [ meta.id, file(meta.read1), file(meta.read2) ] 
        }
    }
    .set { ch_reads }

    // get crams that need to be converted to fastq
    ch_sample_meta
    .map { meta -> 
        if ((params.rundir && meta.i5index && meta.i7index) || meta.fastq_list) {
            if (meta.dragen_path) {
                if (file(meta.dragen_path).isDirectory()) {
                    def cram = file("${meta.dragen_path}/${meta.id}_tumor.cram", checkExists: true)
                    [ meta.id, cram ]
                } else if (file(meta.dragen_path).isFile()) {
                    [ meta.id, file(meta.dragen_path) ]
                }
            }
        }
    }
    .set { ch_cramtofastq }


    MAKE_FASTQLIST ( ch_fastqlists )

    MAKE_FASTQLIST.out.fastqlist
    .splitCsv ( header:true, sep:',')
    .map { id, fastqlist, runinfo -> [ id, parse_fastqlist(fastqlist), runinfo ] }
    .mix(ch_fastqs)
    .set { ch_fastqs }
    ch_versions = ch_versions.mix(MAKE_FASTQLIST.out.versions)

    SPLIT_CRAM(ch_cramtofastq)
    ch_versions = ch_versions.mix(SPLIT_CRAM.out.versions)

    FASTQ_FROM_CRAM(SPLIT_CRAM.out.crams.transpose())      
    
    FASTQ_FROM_CRAM.out.fastqs
    .concat(ch_reads)
    .set { ch_reads }
    ch_versions = ch_versions.mix(FASTQ_FROM_CRAM.out.versions)

    GET_FASTQ_INFO(ch_reads)

    GET_FASTQ_INFO.out.fastqlist
    .splitCsv ( header:true, sep:',')
    .map { id, fastqlist, runinfo -> [ id, parse_fastqlist(fastqlist), runinfo ] }
    .mix(ch_fastqs)
    .set { ch_fastqs }
    ch_versions = ch_versions.mix(GET_FASTQ_INFO.out.versions)

emit:
    fastqlist = ch_fastqs
    versions = ch_versions
}


def parse_fastqlist(LinkedHashMap row) {
       
    def read1 = new File(row.Read1File)
    def read2 = new File(row.Read2File)

    if (!file(read1).exists()) { 
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${read1}"
    }
    if (!file(read2).exists()) { 
        exit 1, "ERROR: Please check input samplesheet -> Read 1 FastQ file does not exist!\n${read2}"
    }

    fastq_meta = [ row.RGID, row.RGLB, row.Lane, file(read1), file(read2) ]

    return fastq_meta

}
