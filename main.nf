#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Default parameter values

OUTDIR = "${params.outdir}"




log.info """\
         E X O M E - R N A S E Q - N F   P I P E L I N E
         ===================================
         NF version   : $nextflow.version
         runName      : $workflow.runName
         username     : $workflow.userName
         configs      : $workflow.configFiles
         cmd line     : $workflow.commandLine
         start time   : $workflow.start
         projectDir   : $workflow.projectDir
         launchDir    : $workflow.launchDir
         workdDir     : $workflow.workDir
         homeDir      : $workflow.homeDir
         samplesheet  : ${params.samplesheet}
         outdir       : ${params.resultsdir}
         """
         .stripIndent()



process PREPARE_SAMPLESHEET {

    input:
    path samplesheet
    val genome_version

    output:
    path("*csv")

    script:
    """

    if [ ${genome_version} == "hg19" ]; then
        python ${workflow.projectDir}/bin/split_samplesheet.py ${samplesheet} .
    elif [ ${genome_version} == "mm39" ]; then
        cp ${samplesheet} mouse_rnaseq.csv
    else
        echo "Error: Unknown genome: ${genome_version}"
        exit 1
    fi
    """
}

include {RNAseq_only} from './workflows/RNAseq_only.nf'
include {RNAseq_multiple_libs} from './workflows/RNAseq_multiple_libs.nf'
include {Exome_only_WF} from './workflows/Exome_only_WF.nf'
include {Tumor_multiple_libs} from './workflows/Tumor_multiple_libs.nf'
include {Tumor_Normal_WF} from './workflows/Tumor_Normal_WF.nf'
include {Tumor_Normal_RNAseq_WF} from './workflows/Tumor_Normal_RNAseq_WF.nf'
include {Tumor_RNAseq_WF} from './workflows/Tumor_RNAseq_WF.nf'
include {Mouse_RNA} from './workflows/Mouse_RNA.nf'


workflow {
prepared_samplesheets = PREPARE_SAMPLESHEET(params.samplesheet,params.genome_v)


    prepared_samplesheets.branch {
        rnaseq: it.name == 'RNAseq.csv'
        exome: it.name == 'Exome.csv'
        multiple_exome: it.name == 'Tumor_lib.csv'
        tumor_rnaseq_normal: it.name == 'Tumor_RNAseq_Normal.csv'
        multiple_rna: it.name == 'RNA_lib.csv'
        tumor_normal: it.name == 'Tumor_Normal.csv'
        tumor_rnaseq: it.name == 'Tumor_RNAseq.csv'
        mouse_rna: it.name == 'mouse_rnaseq.csv'


    }.set { branched_samplesheets }

    branched_samplesheets.rnaseq | RNAseq_only
    branched_samplesheets.tumor_rnaseq | Tumor_RNAseq_WF
    branched_samplesheets.exome | Exome_only_WF
    branched_samplesheets.multiple_exome | Tumor_multiple_libs
    branched_samplesheets.tumor_rnaseq_normal | Tumor_Normal_RNAseq_WF
    branched_samplesheets.multiple_rna | RNAseq_multiple_libs
    branched_samplesheets.tumor_normal | Tumor_Normal_WF
    branched_samplesheets.mouse_rna | Mouse_RNA

}

workflow.onComplete {

    def message = Utils.handleWorkflowCompletion(
        workflow,
        "biowulf_mouse_RNA_slurm",
        workflow.profile == "biowulf_mouse_RNA_slurm" ? "completed.txt" : "successful.txt"

    )
    if (workflow.profile.contains("biowulf")) {
        sendMail(
            to: "${workflow.userName}@mail.nih.gov",
            cc: workflow.profile == "biowulf_mouse_RNA_slurm" ? "" : "wenxi@mail.nih.gov, gangalapudiv2@mail.nih.gov",
            subject: workflow.success ? "khanlab ngs-pipeline execution successful" : "khanlab ngs-pipeline execution failed",
            body: message,
            mimeType: 'text/html'
        )
    }
}
