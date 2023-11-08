process Sequenza_utils {

     tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/sequenza", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(Nbam),path(Nindex),
     path(Tbam),path(Tindex),path(bed),
     path(genome),path(genome_fai),
     path(genome_dict),path(gc50base)

     output:
     tuple val(meta),path("${meta.lib}.seqz_small.gz") , emit: sequenza_bin
     path "versions.yml"             , emit: versions

    stub:
    """
        touch "${meta.lib}.seqz_small.gz"

    """

    script:

    """
    sequenza-utils bam2seqz --normal ${Nbam} --tumor ${Tbam} --fasta ${genome} -gc ${gc50base} -o ${meta.lib}.seqz.gz
    sequenza-utils seqz_binning -w 50 -s ${meta.lib}.seqz.gz | gzip > ${meta.lib}.seqz_small.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sequenza-utils: \$(sequenza-utils --version 2>&1 |grep -E '^This is version'|cut -f4 -d " ")
    END_VERSIONS
    """
}

process Sequenza {

     tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/sequenza", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(seqz_small)

     output:
     tuple val(meta), path("${meta.lib}/${meta.lib}_segments.txt"), emit: segments
     tuple val(meta), path("${meta.lib}/${meta.lib}*pdf"), emit: pdf
     tuple val(meta), path("${meta.lib}/${meta.lib}_mutations.txt"), emit: mutations
     path("${meta.lib}/${meta.lib}_alternative_solutions.txt"), emit: alternate
     tuple val(meta), path("${meta.lib}/${meta.lib}_confints_CP.txt"), emit: confints

    stub:
    """
    touch "${meta.lib}/${meta.lib}_segments.txt"

    """

    script:

    """
    run_sequenza_pipeline.R  --sample ${meta.lib}

    """
}


process Sequenza_annot {
     tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/sequenza", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(segments),path(capture_bed),path(combined_gene_list)

     output:
     tuple val(meta),
     path("${meta.lib}.txt")

    stub:
    """
    touch "${meta.lib}.txt"

    """

    script:

    """
     set +eo pipefail
     sed -e 's/"//g' ${segments} |sed -e 's/chromosome/#chromosome/' |  bedtools intersect -wa -a - -b ${capture_bed} |grep -v NOTFOUND |sed -e 's/___/\t/g'| cut -f 1-4|bedtools expand -c 4 > ${meta.lib}.txt.tmp
     sed -e 's/"//g' ${segments} |sed -e 's/chromosome/#chromosome/' |head -1 > ${meta.lib}.txt.tmp1
     sed -i 's/end.pos\tBf/end.pos\tGene\tBf/g' ${meta.lib}.txt.tmp1
     sed -e 's/"//g' ${segments} |sed -e 's/chromosome/#chromosome/' |intersectBed -a - -b ${meta.lib}.txt.tmp -wb|cut -f 1-4,8-100 >>${meta.lib}.txt.tmp1
     GeneAnnotation.v1.pl ${combined_gene_list} ${meta.lib}.txt.tmp1 3 > ${meta.lib}.txt
     rm -rf ${meta.lib}.txt.tmp1 ${meta.lib}.txt.tmp
    """






}
