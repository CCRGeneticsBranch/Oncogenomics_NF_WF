process fastqc {
        tag { dataset_id }
        publishDir "$params.resultsdir/$dataset_id", mode: 'copy'

        cache false
        input:
        tuple val(dataset_id),
        path(forward),
        path(reverse)

        output:
        tuple val("${dataset_id}"),
        path("fastqc_trim_${dataset_id}")

        container 'nciccbr/ccbr_fastqc_0.11.9:v1.1'

        script:
        """
        mkdir fastqc_trim_${dataset_id}
        fastqc -o fastqc_trim_${dataset_id} -q $forward $reverse
        """
}

