process Fastqc {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc/", mode: 'copy'
    
    input:
    tuple val(meta), path(trim), path(r1fq), path(r2fq)
    
    output:
    tuple val(meta), path("fastqc") 
    
    script:
    def args = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.lib}"
    
    """
    if [ ! -d fastqc ];then mkdir -p fastqc;fi
    fastqc  ${trim[0]} ${trim[1]} $r1fq $r2fq -t $task.cpus -o fastqc
    """
}


process Multiqc {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path("*")
 

    output:
    tuple val(meta),path("multiqc_report.html")
    //path "versions.yml"

    script:
    """
    multiqc . -f

    
    """

}

process Genotyping {
    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(bam),
        path(index),
        path(Sites1000g4genotyping),
        path(genome),
        path(genome_fai),
        path(genome_dict)

    output:
    tuple val(meta),
    path("${meta.lib}.gt"),
    path("${meta.lib}.loh")

    stub:
    """
    touch "${meta.lib}.gt"
    touch "${meta.lib}.loh"
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
     """
    bcftools mpileup -R ${Sites1000g4genotyping} -C50 -Oz -d 1000 -f ${genome} ${bam} | bcftools call --ploidy GRCh37 -mv -Ov -o ${prefix}.samtools.vcf

    vcf2genotype.pl ${prefix}.samtools.vcf > ${prefix}.gt

    vcf2loh.pl ${prefix}.samtools.vcf  > ${prefix}.loh


    """
}


process CircosPlot {
    
    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(gt),
        path(loh)

    output:
    tuple val(meta),
        path("${meta.lib}.circos.png")

    stub:
    """
    touch "${meta.lib}.circos.png"
    """

    script:

     """
     circosLib.R  \$PWD/ ${meta.lib}.circos.png ${meta.lib}
     """
}


process RNAseQC {

    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(rRNA_interval),
        path(transcript_gtf)

    output:
    tuple val(meta),
        path("rnaseqc/report.html")

    stub:
     """
     touch "report.html"
     """

    script:
     """
     java -jar \$RNASEQCJAR -r ${genome} -rRNA ${rRNA_interval} -o rnaseqc -s "${meta.lib}|${bam}|${meta.lib}" -t ${transcript_gtf}
 
     """

}

