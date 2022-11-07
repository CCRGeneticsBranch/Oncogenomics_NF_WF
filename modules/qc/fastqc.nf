process fastqc {
        tag { dataset_id }
        publishDir "$params.resultsdir/$dataset_id", mode: 'copy'

        input:
        tuple val(dataset_id),
        path(forward),
        path(reverse)

        output:
        tuple val("${dataset_id}"),
        path("fastqc_trim_${dataset_id}")


        script:
        """
        mkdir fastqc_trim_${dataset_id}
        fastqc -o fastqc_trim_${dataset_id} -q $forward $reverse
        """
}

