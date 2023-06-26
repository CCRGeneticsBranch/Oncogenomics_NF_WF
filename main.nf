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

include {RNAseq_only} from './workflows/RNAseq_only.nf'
include {RNAseq_multiple_libs} from './workflows/RNAseq_multiple_libs.nf'
include {Tumor_only_WF} from './workflows/Tumor_only_WF.nf'
include {Tumor_multiple_libs} from './workflows/Tumor_multiple_libs.nf'


workflow {
/*   // Check if Tumor_lib.csv is present
    if (fileExists("Tumor_lib.csv")) {
        // Launch Exome workflow
        Tumor_multiple_libs()
    } else {
        // Print a message indicating that Tumor_lib.csv is not present
        println("No Exome workflow to run. Tumor_lib.csv is missing.")
    }
*/
    if (fileExists("Tumor.csv")) {
        // Launch Exome workflow
        Tumor_only_WF()
    } else {
        // Print a message indicating that Tumor_lib.csv is not present
        println("No Exome workflow to run. Tumor_lib.csv is missing.")
    }

    // Check if RNAseq.csv is present

    if (fileExists("RNAseq.csv")) {
        // Launch RNASEQ workflow
        RNAseq_only()
    } else {
        // Print a message indicating that RNAseq.csv is not present
        println("No RNASEQ workflow to run. RNAseq.csv is missing.")
    }


    if (fileExists("RNA_lib.csv")) {
        // Launch RNASEQ workflow
        RNAseq_multiple_libs()
    } else {
        // Print a message indicating that RNAseq.csv is not present
        println("No RNASEQ workflow to run. RNAseq.csv is missing.")
    }


}

def fileExists(String filePath) {
    new File(filePath).exists()
}

