process rsem {
        tag { dataset_id }

        publishDir "$params.resultsdir/$dataset_id", mode: 'move'

        input:
        tuple val(dataset_id),
        path(bam),
        path(genomeIndex)

        output:
        tuple val("${dataset_id}"), path("trim_${dataset_id}.genes.results")
        
        container 'docker://nciccbr/ccbr_rsem_1.3.3:v1.0'

        script:
        """
        rsem-calculate-expression --no-bam-output --paired-end -p ${task.cpus}  --estimate-rspd  --bam $bam ${genomeIndex}/rsem_1.3.2 trim_$dataset_id
        """

}

