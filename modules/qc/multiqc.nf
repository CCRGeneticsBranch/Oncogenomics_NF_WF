process multiqc {
        tag { dataset_id }
//        publishDir "$params.resultsdir/$dataset_id", mode: 'move'
        publishDir "s3://agc-424336837382-us-east-1/nfmvpout/$dataset_id", mode: 'move'
        input:
        tuple val(dataset_id),
        path(qc)
        output:
        path "multiqc_report.html"
//      tuple val("${dataset_id}"), path("trim_${dataset_id}_multiqc_report.html")
        container 'nciccbr/ccbr_multiqc_1.9:v0.0.1'

        script:
        """
        multiqc -m fastqc .

        """

}

