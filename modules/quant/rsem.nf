process rsem {
        tag { dataset_id }

//        publishDir "s3://agc-424336837382-us-east-1/nfmvpout/$dataset_id", mode: 'copy'
        publishDir "$params.resultsdir/$dataset_id", mode: 'copy'

        input:

        tuple val(dataset_id),
        path(T_bam),
        path(G_bam),
        path(genomeIndex)

        output:
        tuple val("${dataset_id}"), path("trim_${dataset_id}.genes.results")
        
        container 'nciccbr/ccbr_rsem_1.3.3:v1.0'

        script:
        """
        rsem-calculate-expression --no-bam-output --paired-end -p ${task.cpus}  --estimate-rspd  --bam $T_bam ${genomeIndex}/rsem_1.3.2 trim_$dataset_id
        """

}

