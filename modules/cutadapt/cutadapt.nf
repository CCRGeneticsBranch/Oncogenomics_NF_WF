process Cutadapt {
     tag "$meta.lib"
//        publishDir "$params.resultsdir/${meta.id}/${meta.casename}/${meta.lib}/Cutadapt", mode: 'copy'
     input:
     tuple val(meta), path(r1fq), path(r2fq)

     output:
     tuple val(meta), path("*.gz") , emit: trim_reads
     path "versions.yml"             , emit: versions

     script:
     def args = task.ext.args   ?: ''
     def prefix   = task.ext.prefix ?: "${meta.lib}"

     """
     cutadapt --pair-filter=any \\
        --nextseq-trim=2 \\
        --trim-n \\
        -n 5 -O 5 \\
        -q 10,10 -m 30:30 \\
        -b file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa \\
        -B file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa \\
        -j $task.cpus \\
        -o ${prefix}_R1.trim.fastq.gz \\
        -p ${prefix}_R2.trim.fastq.gz \\
        $r1fq $r2fq

     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
         Cutadapt: \$(cutadapt --version)
     END_VERSIONS
     """

}
