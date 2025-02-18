
import nextflow.Nextflow

class HandleWorkflowCompletion {

    static String handleWorkflowCompletion(workflow, String profileName, String successFileName) {
        def workDir = new File("${workflow.launchDir}/work")

        // Update permissions
        if (workDir.exists()) {
            println "Updating permissions for work directory: ${workDir}"
            "chmod -R 775 ${workDir}".execute().waitFor()
        } else {
            println "Work directory not found: ${workDir}"
        }

        def successFilePath = "${workflow.launchDir}/${successFileName}"
        def successFile = new File(successFilePath)

        String fullMessage = ""

        if (workflow.success) {
            // Ensure the success file is overwritten
            if (successFile.exists()) {
                successFile.delete()
            }
            successFile.createNewFile()

            if (workflow.profile == profileName) {
                fullMessage = "Mouse RNAseq workflow completed successfully. Results are located at ${workflow.launchDir}"
            } else {
                def htmlFile = new File("${workflow.launchDir}/qc/genotyping.html")
                def htmlContent = htmlFile.exists() ? htmlFile.text : "QC report not found."
                fullMessage = htmlContent
            }
        } else {
            fullMessage = "Workflow completed with errors. Error log is located at ${workflow.launchDir}"
        }

        return fullMessage
    }
}
