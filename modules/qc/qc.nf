process Fastqc {
    tag { dataset_id }

    publishDir "$params.resultsdir/$dataset_id/qc/", mode: 'copy'

    input:
    tuple val(dataset_id),
        path(r1),
        path(r2),
        path(trim_r1),
        path(trim_r2)

    output:
    tuple val("$dataset_id"),
         path("fastqc")

    script:
    """
    if [ ! -d fastqc ];then mkdir -p fastqc;fi
    fastqc $r1 $r2 $trim_r1 $trim_r2 -t $task.cpus -o fastqc
    """
}


process Multiqc {
    tag { dataset_id }
    cache false
    publishDir "${params.resultsdir}/${dataset_id}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
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
    publishDir "${params.resultsdir}/${dataset_id}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        path(bam),
        path(index),
        path(Sites1000g4genotyping),
        path(genome),
        path(genome_fai),
        path(genome_dict)

    output:
    tuple val("${dataset_id}"),
    path("${dataset_id}.star.gt"),
    path("${dataset_id}.star.loh")

    stub:
    """
    touch "${dataset_id}.star.gt"
    touch "${dataset_id}.star.loh"
    """

    shell:
     '''

    bcftools mpileup -R !{Sites1000g4genotyping} -C50 -Oz  -f !{genome} !{bam} | bcftools call --ploidy GRCh37 -mv -Ov -o !{dataset_id}.star.samtools.vcf

    vcf2genotype.pl !{dataset_id}.star.samtools.vcf > !{dataset_id}.star.gt

    vcf2loh.pl !{dataset_id}.star.samtools.vcf  > !{dataset_id}.star.loh


    '''
}

process CircosPlot {
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        path(gt),
        path(loh)

    output:
    tuple val("${dataset_id}"),
        path("${dataset_id}.star.circos.png")    

    stub:
    """
    touch "${dataset_id}.star.circos.png"
    """

    shell:
     '''
     circosLib.R  $PWD/ !{dataset_id}.star.circos.png !{dataset_id}
     '''
}


process RNAseQC {

    publishDir "${params.resultsdir}/${dataset_id}/qc", mode: "${params.publishDirMode}"

    tag { dataset_id }

    input:
    tuple val(dataset_id),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(rRNA_interval),
        path(transcript_gtf)

    output:
    tuple val("${dataset_id}"),
        path("rnaseqc/report.html")

    stub:
     """
     touch "report.html"
     """

    shell:
     '''
     java -jar $RNASEQCJAR -r !{genome} -rRNA !{rRNA_interval} -o rnaseqc -s "!{dataset_id}|!{bam}|!{dataset_id}" -t !{transcript_gtf}
 
     '''

}
