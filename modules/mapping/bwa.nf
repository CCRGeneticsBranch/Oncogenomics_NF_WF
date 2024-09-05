process BWA {
    tag "$meta.lib"
    //publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}", mode: "${params.publishDirMode}"

    input:
    tuple val(meta), path(r1fq),path(r2fq),path(bwa_genomeindex)

    output:
    tuple val(meta), path("${meta.lib}.bam"), emit: bwa_bam
    tuple val(meta), path("${meta.lib}.bam.bai"), emit: bwa_bai
    path "versions.yml"             , emit: versions

    stub:
    """
    touch "${meta.lib}.bwa.bam"
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"

    """
    bwa mem -M -t ${task.cpus} -R '@RG\\tID:${prefix}\\tSM:${prefix}\\tLB:${prefix}\\tPL:illumina' ${bwa_genomeindex}/hg19.fa ${r1fq} ${r2fq} | samtools view -Sbh - | samtools sort -m 30000000000 -o ${prefix}.bam
    samtools index ${prefix}.bam
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 |grep -E '^Version'|sed 's/.*Version: //')
        samtools: \$(samtools --version |grep -E '^samtools'|sed 's/.*samtools //')
    END_VERSIONS
    """


}
