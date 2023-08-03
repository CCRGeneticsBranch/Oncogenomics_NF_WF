process Sequenza_utils {

     tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/sequenza", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(Tbam),path(Tindex),path(Tbed)
     tuple val(meta2),path(Nbam),path(Nindex)
     path genome
     path genome_fai
     path genome_dict
     path gc50base
     
     output:
     tuple val(meta),
     path("${meta.lib}.seqz_small.gz")

    stub:
    """   
        touch "${meta.lib}.seqz_small.gz"

    """

    script:

    """
    sequenza-utils bam2seqz --normal ${Nbam} --tumor ${Tbam} --fasta ${genome} -gc ${gc50base} -o ${meta.lib}.seqz.gz
    sequenza-utils seqz_binning -w 50 -s ${meta.lib}.seqz.gz | gzip > ${meta.lib}.seqz_small.gz

    """
}

process Sequenza {

     tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/sequenza", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(seqz_small)
     path sequenza_Rscript
     
     output:
     tuple val(meta),
     path("${meta.lib}_CN_bars.pdf")

    stub:
    """   
    touch "${meta.lib}_CN_bars.pdf"

    """

    script:

    """
     ${sequenza_Rscript} --sample ${meta.lib}
    """
}