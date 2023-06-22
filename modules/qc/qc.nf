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
    fastqc --extract ${trim[0]} ${trim[1]} $r1fq $r2fq -t $task.cpus -o fastqc
    """
}

process Fastq_screen {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc/fastq_screen/", mode: 'copy'
    
    input:
    tuple val(meta),path(trim),path(fastq_screen_config),path(fqs_human)
    
    stub:
    """
    touch "${meta.lib}_R1_screen.html"
    touch "${meta.lib}_R1_screen.html"
    """

    output:
    tuple val(meta), 
    path("${meta.lib}_R1_screen.html"),
    path("${meta.lib}_R2_screen.html")

    
    script:
    def args = task.ext.args   ?: ''
    def prefix   = task.ext.prefix ?: "${meta.lib}"
    
    """
    #--subset 1000000  add this parameter after testing.
    if [ ! -d fastq_screen ];then mkdir -p fastq_screen;fi
    fastq_screen --conf ${fastq_screen_config}  --aligner bowtie2 --force ${trim[0]} ${trim[1]}  --outdir fastq_screen
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
    tuple val(meta),path("${meta.lib}.gt") , emit: gt
    tuple val(meta),path("${meta.lib}.loh"), emit: loh

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


process CircosPlot_lib {
    
    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
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

