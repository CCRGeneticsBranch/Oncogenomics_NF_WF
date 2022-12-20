process Rsem {
        tag { dataset_id }

        publishDir "$params.resultsdir/$dataset_id/RSEM", mode: 'copy'

        input:

        tuple val(dataset_id),
        path(T_bam),
        path(G_bam),
        path(G_bai),
        path(chimeric_junction),
        path(genomeIndex)

        output:
        tuple val("${dataset_id}"),
		path("${dataset_id}.genes.results")
        

        script:
        """
        rsem-calculate-expression --no-bam-output --paired-end -p ${task.cpus} --estimate-rspd --bam $T_bam ${genomeIndex}/rsem_1.3.2 $dataset_id
        """

}

