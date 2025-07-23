
import nextflow.Nextflow

class HandleWorkflowCompletion {

    static String handleWorkflowCompletion(workflow, String genome, String successFileName, String platform, String patientId, String casename, String resultsdir) {
        if (platform == "biowulf") {
            println "PatientID: ${patientId}"
            def workDir = new File("${workflow.launchDir}/work")

        // Update permissions
            if (workDir.exists()) {
                println "Updating permissions for work directory: ${workDir}"
                "chmod -R 775 ${workDir}".execute().waitFor()
            } else {
                println "Work directory not found: ${workDir}"
            }
        } else {
            println "Platform does not match. Skipping permission update."
        }

        //def successFilePath = "${workflow.launchDir}/${successFileName}"
        //def successFilePath = "${resultsdir}/${patientId}/${casename}/${successFileName}"
        //def successFile = new File(successFilePath)

        String fullMessage = ""

        if (workflow.success) {
            // Ensure the success file is overwritten
            //if (successFile.exists()) {
            //    successFile.delete()
            //}
            //successFile.createNewFile()

            if (genome == "mm39") {
                fullMessage = "Mouse RNAseq workflow completed successfully. Results are located at ${resultsdir}/${patientId}/${casename}"
            } else {
                def htmlFile = new File("${resultsdir}/${patientId}/${casename}/qc/genotyping.html")
                def htmlContent = htmlFile.exists() ? htmlFile.text : "QC report not found for ${resultsdir}/${patientId}/${casename}"
                fullMessage = htmlContent
            }
        } else {
            fullMessage = "Workflow completed with errors. Error log is located at ${resultsdir}/${patientId}/${casename}"
        }

        return fullMessage
    }
}
