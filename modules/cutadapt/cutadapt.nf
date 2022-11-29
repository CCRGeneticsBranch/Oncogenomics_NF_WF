process cutadapt {
        tag { dataset_id }
//        publishDir 's3://agc-424336837382-us-east-1/nfmvpout', mode: 'copy'
        publishDir "$params.resultsdir/$dataset_id", mode: 'copy'

        input:
        tuple val(dataset_id),
        path(r1fq),
        path(r2fq)

        output:
        tuple val("${dataset_id}"),
        path("${dataset_id}_R1.trim.fastq.gz"),
        path("${dataset_id}_R2.trim.fastq.gz")


        script:
        """
	cutadapt --pair-filter=any \\
	--nextseq-trim=2 \\
	--trim-n \\
	-n 5 -O 5 \\
	-q 10,10 -m 30:30 \\
	-b file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa \\
	-B file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa \\
	-j $task.cpus \\
        -o ${dataset_id}_R1.trim.fastq.gz \\
        -p ${dataset_id}_R2.trim.fastq.gz \\
	$r1fq $r2fq      
        """

}
