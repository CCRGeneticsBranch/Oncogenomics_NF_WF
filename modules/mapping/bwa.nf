process BWA {
    tag "$meta.lib"
    //publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}", mode: "${params.publishDirMode}"

    input:
    tuple val(meta), path(r1fq),path(r2fq),path(bwa_genomeindex)

    output:
    tuple val(meta), 
    path("${meta.lib}.bam"), 
    path("${meta.lib}.bam.bai")

    stub:
    """
    touch "${meta.lib}.bwa.bam"
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"

    """
    bwa mem -M -t ${task.cpus} -R '@RG\\tID:${prefix}\\tSM:${prefix}\\tLB:${prefix}\\tPL:illumina' ${bwa_genomeindex}/hg19.fa ${r1fq} ${r2fq} | samtools view -Sbh - | samtools sort -m 30000000000 -o ${prefix}.bam 
    samtools index ${prefix}.bam
    """


}