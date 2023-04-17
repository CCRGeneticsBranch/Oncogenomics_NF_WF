process Fastqc {
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/${library}/qc/", mode: 'copy'

    input:
    tuple val(dataset_id),
        val(library),
        path(r1),
        path(r2),
        path(trim_r1),
        path(trim_r2)

    output:
    tuple val("$dataset_id"),
         val("$library"),
         path("fastqc")

    script:
    """
    if [ ! -d fastqc ];then mkdir -p fastqc;fi
    fastqc $r1 $r2 $trim_r1 $trim_r2 -t $task.cpus -o fastqc
    """
}


process Multiqc {
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        val(library),
        path("*")
 

    output:
    path("multiqc_report.html")

    script:
    """
    multiqc . -f
    """

}

process Genotyping {
    tag { dataset_id }
    publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/${library}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        val(library),
        path(bam),
        path(index),
        path(Sites1000g4genotyping),
        path(genome),
        path(genome_fai),
        path(genome_dict)

    output:
    tuple val("${dataset_id}"),
    val("${library}"),
    path("${library}.star.gt"),
    path("${library}.star.loh")

    stub:
    """
    touch "${library}.star.gt"
    touch "${library}.star.loh"
    """

    shell:
     '''

    bcftools mpileup -R !{Sites1000g4genotyping} -C50 -Oz -d 1000 -f !{genome} !{bam} | bcftools call --ploidy GRCh37 -mv -Ov -o !{library}.star.samtools.vcf

    vcf2genotype.pl !{library}.star.samtools.vcf > !{library}.star.gt

    vcf2loh.pl !{library}.star.samtools.vcf  > !{library}.star.loh


    '''
}


process CircosPlot {
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        val(library),
        path(gt),
        path(loh)

    output:
    tuple val("${dataset_id}"),
        val("${library}"),
        path("${dataset_id}.circos.png")

    stub:
    """
    touch "${dataset_id}.circos.png"
    """

    shell:
     '''
     circosLib.R  $PWD/ !{dataset_id}.circos.png !{library}
     '''
}


process RNAseQC {

    publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/${library}/qc", mode: "${params.publishDirMode}"

    tag { dataset_id }

    input:
    tuple val(dataset_id),
        val(library),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(rRNA_interval),
        path(transcript_gtf)

    output:
    tuple val("${dataset_id}"),
        val("${library}"),
        path("rnaseqc/report.html")

    stub:
     """
     touch "report.html"
     """

    shell:
     '''
     java -jar $RNASEQCJAR -r !{genome} -rRNA !{rRNA_interval} -o rnaseqc -s "!{library}|!{bam}|!{library}" -t !{transcript_gtf}
 
     '''

}
