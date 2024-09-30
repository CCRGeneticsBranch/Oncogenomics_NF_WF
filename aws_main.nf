nextflow.enable.dsl=2

log.info """\
         E X O M E - R N A S E Q - N F   P I P E L I N E
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


process Check_samplesheet {

	publishDir "${params.resultsdir}"

	input:
	path(input_csv)

	output:
	path("*csv")

	script:
	"""
	${workflow.projectDir}/bin/split_samplesheet.py ${input_csv} .
	"""

}

include {RNAseq_only} from './workflows/RNAseq_only.nf'

workflow {
    Check_samplesheet(params.samplesheet)
	basename = Check_samplesheet.out
            .map{ file -> file.baseName}
    basename.view()


    if (basename.contains("RNAseq")) {
        RNAseq_only()
    } else {
        println("No workflow to run. Required file is missing.")
    }


}
