process BATCH_QC {
  container "registry.gsc.wustl.edu/mgi-cle/chromoseq:v1.1"
  memory '4 GB'

  input:
  tuple val(meta), val (orderby)
  val qcOut = ${params.outdir} + "/QC_info.txt"

  output:
  stdout

  script: 
    """
    if [ -n "$(/bin/ls -d ${params.outdir}/TWDY-*)" ]; then
        /bin/chmod -R 777 ${params.outdir}/TWDY-*
    fi

    /usr/bin/perl /usr/local/bin/QC_info.pl "${params.outdir}/*/*.chromoseq.txt" > ${qcOut}
    """
}