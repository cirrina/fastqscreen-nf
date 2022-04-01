#!/bin/bash -l


#######################################
# driver for  fastqscreen standalone pipe
#######################################

## INPUT: CTG IEM Sample Sheet
## INPUT: FASTQ files (default fastq parth is output from bcl2fastq analysis)


# scripts_root='/Users/david/scripts' #

# script execution dir. If generate workfolder & copy scripts, this WILL be checked & warn against the script version rovided in samplesheet
# if initiating script in project workfolder, the script here is used & version will not be checked
script_exec_dir=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd ) # Where is this script located...
script_exec_dir=$(cd ${script_exec_dir} && pwd -P) # needed in case script exec dir is a symlink
## Root directories (based on ourr set folder naming conventions on lsens)

scripts_root="/projects/fs1/tools/bin/"
project_root='/projects/fs1/shared/ctg-projects'
delivery_root='/projects/fs1/shared/ctg-delivery' ## should be added pipelineProfile/ProjectID
ctg_qc_root='/projects/fs1/shared/ctg-qc/' ## should be added pipelineProfile/ProjectID

## Script & config names
nf_script="fastqscreen-main.nf"
pipelineName="fastqscreen-nf"


#######################
#  == Initiation ==
#######################
# fastq_custom=true ## set to false if no fastq_input_dir is supplied (will defalt to {project_})
prime_projectfolder_mode='false' ## -p true,  if to generate project folder and configs, but not initiate the nextgflow
resume='false' ## used to resume nextflow pipeline, Need to be executed within project workfolder.
exec_dir=$(pwd)
fastq_input_dir=$(pwd)
projectid=${PWD##*/}
samplesheet=""
echo ${projectid}

# usage message
usage() {
    echo ""
    echo " fastqscreen-driver  \
[ -h help ] "  1>&2
    echo " help                       -h : print help message"

}

exit_abnormal() {
  usage
  exit 1
}

#######################
#  == Initiation ==
#######################
# fastq_custom=true ## set to false if no fastq_input_dir is supplied (will defalt to {project_})
prime_projectfolder_mode='false' ## -p true,  if to generate project folder and configs, but not initiate the nextgflow
resume='false' ## used to resume nextflow pipeline, Need to be executed within project workfolder.
exec_dir=$(pwd)
fastq_input_dir=$(pwd)
projectid=${PWD##*/}
samplesheet=""
echo ${projectid}

# usage message
usage() {
    echo ""
    echo " fastqscreen-driver  \
[ -h help ] "  1>&2
    echo " help                       -h : print help message"

}

exit_abnormal() {
  usage
  exit 1
}



################################################
# == 1 ==   Check input arguments
################################################

while getopts "sh:" opt; do
    case $opt in
      s) samplesheet=$OPTARG
         ;;
      h) exit_abnormal
	;;
      \?) echo echo ""; echo " Error ctg-rnaseq : ";"Invalid option -$OPTARG" >&2
        exit_abnormal ;;
      :) echo ""; echo "Error:"; echo " -${OPTARG} requires an argument!"echo ""; echo ""
                 exit 1 ;;
    esac
done



shift "$(( OPTIND -1 ))"


## if no (or wrongly specified) samplesheet provided, then generate samplesheet from all available fastq files in dir
if [[ ! -f ${samplesheet} ]]; then
  echo "No samplesheet specified, or Sample Sheet does not exist (in current dir)"

  samplesheet="fastq_samplesheet.csv"
  echo "... Generating Samplesheet from available fastq files: ${samplesheet}"
  echo "[Header]" > ${samplesheet}

  echo "ProjectID",${projectid} >> ${samplesheet}
 echo "PairedReads,false" >> ${samplesheet}
 echo "[Data]" >> ${samplesheet}
 echo "Sample_ID,fastq_1,fastq_2" >> ${samplesheet}

 for file in *fastq.gz
 do
   name=$(echo $file | sed "s/.fastq.gz//g"\ )
   echo $name,$file >> ${samplesheet}
 done
fi


paired_global=$(awk -F, '$1 == "PairedReads"' ${samplesheet} | awk -F, '{print $2}')
echo "  ... paired_global, from samplesheet [Header]: ${paired_global}"
if [[ ${paired_global} != "true" ]] && [[ ${paired_global} != "false"  ]]; then
 echo ""; echo ""; echo " ---------------------------  "; echo " ERROR "; echo " ---------------------------  ";echo ""
 echo " paired_global is not properly supplied in samplesheet."; echo " Must be supplied as 'PairedReads' 'true' OR 'false' within the [Header] section of sample sheet"; echo""; echo ""
 exit 1
fi


# set project work dir
project_dir="${project_root}/fastqscreen/${projectid}"
delivery_dir="${delivery_dir}/fastqscreen/${projectid}"
ctg_qc_dir="${ctg_qc_root}/fastqscreen/${projectid}"

nf_config_project="${project_dir}/nextflow.config.params.${projectid}"


 ## Warnings & Prompt.
 echo ""
 echo "  ... ... scripts_dir:  ${script_exec_dir}"
 echo "  ... ... project workfolder:  ${project_dir} "
 echo ""



 ## Create project directory
 ## ------------------------
 mkdir -p ${project_dir}
 ## Copy scripts from version specific pipeline scripts dir
 if [[ ${script_exec_dir} != ${exec_dir} ]]; then
   cp -r ${script_exec_dir}/* ${project_dir}/ # copy all scripts to workfolder. Will overwrite netflow.config
 fi

 ## Copy samplesheet to project workfolder
 cp ${samplesheet} ${project_dir}
 cd ${project_dir}



   # Create nextflow configuration file -- nextflow.params_${ProjectID} --
   #######################################################################


   echo ""
   echo " ... Writing nextflow parameters to project-specific config: ${nf_config_project}"

   ## Write nextflow params to file
  echo ""  > ${nf_config_project}
  echo "//  nextflow configuration file"                              >> ${nf_config_project}
  echo "//  Project:  ${projectid}"                                   >> ${nf_config_project}
  echo ""                                                             >> ${nf_config_project}
  echo "//  project specific parameters"                              >> ${nf_config_project}
  echo "//  will override params in 'nextflow.config' "               >> ${nf_config_project}
  echo "//"                                                           >> ${nf_config_project}
  echo ""                                                             >> ${nf_config_project}
  echo " params {"                                                      >> ${nf_config_project}
  echo ""                                                               >> ${nf_config_project}
  echo "  // Pipeline                                                 " >> ${nf_config_project}
  echo "  pipelineName       =  '${pipelineName}'                      " >> ${nf_config_project}
  echo "  exec_dir    =  '${exec_dir}'                   " >> ${nf_config_project}
  echo "  script_execution_dir  =  '${script_exec_dir}'               " >> ${nf_config_project}
  echo ""                                                               >> ${nf_config_project}
  echo "  // Project and run directories                               " >> ${nf_config_project}
  echo "  projectid          =  '${projectid}'                       " >> ${nf_config_project}
  echo "  project_dir        =  '${project_dir}'                      " >> ${nf_config_project}
  echo "  delivery_dir       =  '${delivery_dir}'                      " >> ${nf_config_project}
  echo "  ctg_qc_dir         =  '${ctg_qc_dir}'                   " >> ${nf_config_project}
  echo "  fastq_input_dir    =  '${fastq_input_dir}'                        " >> ${nf_config_project}
  echo ""                                                               >> ${nf_config_project}
  echo "  paired_global      =   ${paired_global}                     " >> ${nf_config_project}
  echo ""                                                               >> ${nf_config_project}
  echo "  //  samplesheets                                            " >> ${nf_config_project}
  echo "  samplesheet        =  '${samplesheet}'               " >> ${nf_config_project}
  echo ""                                                               >> ${nf_config_project}
  echo "  // root directories                                         " >> ${nf_config_project}
  echo "  project_root       =  '${project_root}'                    " >> ${nf_config_project}
  echo "  delivery_root      =  '${delivery_root}'                   " >> ${nf_config_project}
  echo "  ctg_qc_root      =  '${ctg_qc_root}'                   " >> ${nf_config_project}
  echo ""                                                               >> ${nf_config_project}
  echo " }"                                                             >> ${nf_config_project}
  echo ""                                                               >> ${nf_config_project}

  #  Priming of project folder complete
  ## --------------------------------
  echo "";echo ""
  echo "  Project primed"
  echo "  ... Project dir        :  ${project_dir}"
  echo "  ... Nextflow config    :  ${nf_config_project}"
  echo "  ... samplesheet        :  ${samplesheet}"
  echo "  ... Scripts dir        :  ${script_exec_dir}"
  echo ""

  chmod -R 775 ${project_dir}
  cd ${project_dir}
