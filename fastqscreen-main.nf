
/* ===============================================================
  *      PERFORM FASTQC, FASTQSCREEN AND MULTIQC ON FASTQ FILES
  =============================================================== */



/* ===============================================================
  *      PARAMS FROM CONFIGS
  =============================================================== */

//  re-assign some params from nextflow.configs
// ----------------------------------------------------
//  project and run folders
projectid           =  params.projectid    // ctg project id. e.g. 2021_024
project_dir         =  params.project_dir   // .../shared/ctg-rnaseq/uroscan/2021_024 // NOT to be confused with project_dir
samplesheet         =  params.samplesheet           // name of simple sample sheet used for pipeline. Must includes file paths to fastq and bamsm, as well as species etc.

fastq_input_dir     =  params.fastq_input_dir            // subdirectory fastq files are located. For a default run the output location from blc2fastq. fastq-files will be read according to sample sheet. Defaults to <bcl2fastq_dir>/<projectid>
delivery_dir        =  params.delivery_dir
file(delivery_dir).mkdir() // main nexttlow work dir for output of analyses. Subdirectory of the project foilder. Files and folders will be moved and copiued from this folder upon pipeline  completion.

// ctg_qc_dir          =  params.ctg_qc_dir
// file(ctg_qc_dir).mkdir() // main nexttlow work dir for output of analyses. Subdirectory of the project foilder. Files and folders will be moved and copiued from this folder upon pipeline  completion.
// runfolder = params.runfolder
// runfolder_path      =  params.runfolder_path

/* ===============================================================
  *      DEFINE PATHS FROM INPUT PARAMS
  =============================================================== */


// Process & Module specific paths (Delivered to customer)
samplesheetsdir = delivery_dir+'/samplesheets'
deliveryscripts = delivery_dir+'/scripts'

// QC dirs
qcdir = delivery_dir+'/qc'
multiqcdir  =  qcdir+'/multiqc'
fastqscreendir = qcdir+'/fastqscreen'




/* ===============================================================
  *       create output and logdirs
  =============================================================== */
// log file for nextflow .onComplete
logfile   =  file( project_dir + '/' + 'log.nextflow.complete' )



/* ===============================================================
  *       CHECKS FILES AND PARAMS
  =============================================================== */


//  Check paramters
// -----------------------------
if (projectid         == '') {exit 1, "You must define a project_id in the nextflow.config"}
if (samplesheet   == '') {exit 1, "You must define a sample sheet path in the nextflow.config"}


// Check if files and directories exist
checkPathParamList = [
  project_dir,
  samplesheet
]
for (param in checkPathParamList) {
    if (param) {
	     file(param, checkIfExists: true)
    }
}




/* ===============================================================
  *       MESSAGES
  =============================================================== */

// Define messages to print and for logfiles
def msg_startup = """\

    Workflow execution parameters
    ---------------------------------
    project id              :  ${projectid}
    project work dir        :  ${project_dir}
    project delivery dir    :  ${delivery_dir}
    ctg qc dir              :  ${ctg_qc_dir}
    nextflow work dir       :  ${workDir}
    samplesheet             :  ${samplesheet}
   """
       .stripIndent()
println( msg_startup )


workflow.onComplete {

  def msg_completed = """\
  	Pipeline execution summary
  	---------------------------
  	Completed at : ${workflow.complete}
  	Duration     : ${workflow.duration}
  	Success      : ${workflow.success}
  	scriptFile   : ${workflow.scriptFile}
    exit status  : ${workflow.exitStatus}
  	errorMessage : ${workflow.errorMessage}
  	errorReport  :
  	"""
  	.stripIndent()
  def error = """\
		${workflow.errorReport}
	   """
  logfile.text = msg_startup.stripIndent()
  logfile.append( msg_completed.stripIndent() )
  logfile.append( error )

  println( msg_completed )
}





/* ===============================================================
  *    PROCESS SAMPLE SHEET & DEFINE CHANNELS
  =============================================================== */
// Read and process sample sheet. Save to SampleSheet-nexflow.csv
// samplesheet to be parsed as input channel (take everything below [Data] section).
sheet = file(params.samplesheet)
all_lines = sheet.readLines()
write_row = false // if next lines has sample info
sheet_nf = file("${project_dir}/SampleSheet-nexflow.csv")
sheet_nf.text=""

for ( line in all_lines ) {
  if ( write_row ) {
    sheet_nf.append(line + "\n")
  }
  if (line.contains("[Data]")) {
    write_row = true
  }
}


/* ===============================================================
  *     Define Channels based from SampleSheet
  =============================================================== */
Channel
  .fromPath(sheet_nf)
  .splitCsv(header:true)
  .map { row -> tuple( row.Sample_ID, row.fastq_1, row.fastq_2 ) }
  .tap{ infoall }
  .set { fastq_ch }

println " > Samples to process: "
println "[Sample_ID,fastq1,fastq2 ]"
infoall.subscribe { println "Info: $it" }




/* ===============================================================
  *      FASTQC
  =============================================================== */

process fastqc {
  tag  { params.run_fastqc  ? "${projectid}__${sid}" : "blank_run"  }
  cpus { params.run_fastqc  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_fastqc  ? params.mem_standard : params.mem_min  }

  input:
  set sid, read1, read2 from fastqc_ch  // from move_fastqc_ch

  output:
  val "x" into fastqc_complete_ch
  set sid, read1, read2, species into fastqscreen_ch

  script:
  if ( params.paired_global && params.run_fastqc)
    """
    mkdir -p ${fastqcdir}
    echo "running fastqc in paired reads mode"
    fastqc ${fastq_input_dir}/${read1} ${fastq_input_dir}/${read2}  --outdir ${fastqcdir}
    """
  else if ( !params.paired_global && params.run_fastqc)
    """
    mkdir -p ${fastqcdir}
    echo "running fastqc in non paired reads mode "
    fastqc ${fastq_input_dir}/${read1}  --outdir ${fastqcdir}
    """
  else
    """
    echo "run_fastqc skipped"
    """
}



/* ===============================================================
  *      FASTQSCREEN
  =============================================================== */

process fastqscreen {

    tag  { params.run_fastqscreen  ? "${projectid}__${sid}" : "blank_run"  }
    cpus { params.run_fastqscreen  ? params.cpu_standard : params.cpu_min  }
    memory { params.run_fastqscreen  ?  params.mem_standard : params.mem_min  }


    input:
    val x from fastqc_complete_ch.collect()
    set sid, read1, read2 from fastqscreen_ch //

    output:
    val "x" into fastqscreen_complete_ch

    script:
    if ( params.paired_global ){
        fqsfiles = "${fastq_input_dir}/${read1} ${fastq_input_dir}/${read2}" }
    else{
        fqsfiles = "${fastq_input_dir}/${read1}" }

    if ( params.run_fastqscreen)
      """
      mkdir -p ${fastqscreendir}
      fastq_screen \\
          --conf ${params.fastqscreen_config} \\
          --subset 500000 \\
          --outdir ${fastqscreendir} \\
          ${fqsfiles}
      """
    else
      """
      echo "run_fastqscreen skipped"
      """
}


process multiqc {
  tag  { params.run_multiqc  ? "${projectid}" : "blank_run"  }
  cpus { params.run_multiqc  ? params.cpu_standard : params.cpu_min  }
  memory { params.run_multiqc  ?  params.mem_standard : params.mem_min  }

  input:
  val x from fastqscreen_complete_ch.collect()

  output:
  val "x" into multiqc_complete_ch

  script:
  mqcreport = multiqcdir + '/' + projectid + '_multiqc_report'

  if ( params.run_multiqc )
    """
    ## use -f flag to overwrite if multiqc is already present from failed run.
    cd ${delivery_dir}
    multiqc -n ${mqcreport} \\
      --interactive \\
      -f \\
      -o ${multiqcdir} . ${fastq_input_dir}

    """
  else
    """
    echo "run_multiqc skipped"
    """
}
