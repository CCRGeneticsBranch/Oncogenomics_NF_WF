// Using DSL-2
nextflow.enable.dsl=2
import groovy.json.JsonSlurper

log.info """\
         R N A S E Q - N F   P I P E L I N E  
         ===================================
         NF version   : $nextflow.version
         runName      : $workflow.runName
         username     : $workflow.userName
         configs      : $workflow.configFiles
         profile      : $workflow.profile
         cmd line     : $workflow.commandLine
         start time   : $workflow.start
         projectDir   : $workflow.projectDir
         launchDir    : $workflow.launchDir
         workdDir     : $workflow.workDir
         homeDir      : $workflow.homeDir
         casename     : ${params.casename}
         """
         .stripIndent()

//import workflows
include { RNASEQ } from './workflows/rnaseq'

workflow {
    RNASEQ ()
}

