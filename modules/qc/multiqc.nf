process multiqc {
        tag { dataset_id }
        publishDir "$params.resultsdir/$dataset_id", mode: 'copy'
//        publishDir "s3://agc-424336837382-us-east-1/nfmvpout/$dataset_id", mode: 'copy'
        
        input:
        tuple val(dataset_id),
        path(qc)
        output:
        path "multiqc_report.html"

        script:
        """
        multiqc -m fastqc .

        """

}

