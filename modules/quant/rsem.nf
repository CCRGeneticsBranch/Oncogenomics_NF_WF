process rsem {
        tag { dataset_id }

//        publishDir "s3://agc-424336837382-us-east-1/nfmvpout/$dataset_id", mode: 'copy'
        publishDir "$params.s3_bucket/nfmvpout/test/$dataset_id", mode: 'copy'
        cache false
        input:
        tuple val(dataset_id),
        path(bam),
        path(genomeIndex)

        output:
        tuple val("${dataset_id}"), path("trim_${dataset_id}.genes.results")
        
        container 'nciccbr/ccbr_rsem_1.3.3:v1.0'

        script:
        """
        rsem-calculate-expression --no-bam-output --paired-end -p ${task.cpus}  --estimate-rspd  --bam $bam ${genomeIndex}/rsem_1.3.2 trim_$dataset_id
        """

}

