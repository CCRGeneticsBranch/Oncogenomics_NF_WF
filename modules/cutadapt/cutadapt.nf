process cutadapt {
        tag { dataset_id }
//        publishDir 's3://agc-424336837382-us-east-1/nfmvpout', mode: 'copy'
        publishDir "$params.s3_bucket/nfmvpout/$dataset_id", mode: 'copy'

        input:
        tuple val(dataset_id),
        path(forward),
        path(reverse)

        output:
        tuple val("${dataset_id}"),
        path("trim_${dataset_id}_R1.fastq"),
        path("trim_${dataset_id}_R2.fastq")

        container 'nciccbr/ncigb_cutadapt_v1.18:latest'

        script:
        """
        cutadapt  -o trim_${dataset_id}_R1.fastq -p trim_${dataset_id}_R2.fastq $forward $reverse
        """

}
