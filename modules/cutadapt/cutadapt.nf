process Cutadapt {
        tag { dataset_id }
        publishDir "$params.resultsdir/$dataset_id/${params.casename}/$library/Cutadapt", mode: 'copy'

        input:

        tuple val(dataset_id),
            val(library),
            path(r1fq),
            path(r2fq)

        output:
        tuple val("${dataset_id}"),
            val("${library}"),
            path("${library}_R1.trim.fastq.gz"),
            path("${library}_R2.trim.fastq.gz")


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
        -o ${library}_R1.trim.fastq.gz \\
        -p ${library}_R2.trim.fastq.gz \\
	$r1fq $r2fq      
        """

}
