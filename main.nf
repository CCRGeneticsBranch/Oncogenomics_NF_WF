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
         """
         .stripIndent()



//import workflows

include {RNAseq_only} from './workflows/RNAseq_only.nf'
include {RNAseq_multiple_libs} from './workflows/RNAseq_multiple_libs.nf'
include {Exome_only_WF} from './workflows/Exome_only_WF.nf'
include {Tumor_multiple_libs} from './workflows/Tumor_multiple_libs.nf'
include {Tumor_Normal_WF} from './workflows/Tumor_Normal_WF.nf'
include {Tumor_Normal_RNAseq_WF} from './workflows/Tumor_Normal_RNAseq_WF.nf'
include {Tumor_RNAseq_WF} from './workflows/Tumor_RNAseq_WF.nf'
include {Mouse_RNA} from './workflows/Mouse_RNA.nf'


// Launch workflow by checking the samplesheet availability
workflow {


if (fileExists("Exome.csv")) {
    Exome_only_WF()
} else if (fileExists("Tumor_lib.csv")) {
    Tumor_multiple_libs()
} else if (fileExists("Tumor_RNAseq_Normal.csv")) {
    Tumor_Normal_RNAseq_WF()
} else if (fileExists("RNAseq.csv")) {
    RNAseq_only()
} else if (fileExists("RNA_lib.csv")) {
    RNAseq_multiple_libs()
} else if (fileExists("Tumor_Normal.csv")) {
    Tumor_Normal_WF()
} else if (fileExists("Tumor_RNAseq.csv")) {
    Tumor_RNAseq_WF()
} else if (fileExists("mouse_rnaseq.csv")) {
    Mouse_RNA()
} else {
    println("No workflow to run. Required file is missing.")
}

}


workflow.onComplete {


    if (workflow.success) {

        def htmlFile = new File("${workflow.launchDir}/qc/genotyping.html")
        def htmlContent = htmlFile.text
        def fullMessage = "${htmlContent}"


        //sendMail(to: "${workflow.userName}@mail.nih.gov" , subject: 'khanlab ngs-pipeline execution', body: fullMessage, mimeType: 'text/html')

        sendMail(to: "${workflow.userName}@mail.nih.gov" , cc: 'khanjav@mail.nih.gov, weij@mail.nih.gov, wenxi@mail.nih.gov, gangalapudiv2@mail.nih.gov,' , subject: 'khanlab ngs-pipeline execution', body: fullMessage, mimeType: 'text/html')

    } else {
        println "Workflow completed with errors. No success email sent."
    }

}


def fileExists(String filePath) {
    new File(filePath).exists()
}
