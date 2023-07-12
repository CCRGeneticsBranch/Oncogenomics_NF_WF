process Actionable_RNAseq {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/Actionable", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
      path(rnaseq),
      path(annotation_coding_rare),
      path(combined_gene_list),
      path(somatic_actionable_sites)

     output:
     tuple val(meta),
        path("${meta.id}.rnaseq.actionable.txt")

     stub:
     """
     touch "${meta.id}.rnaseq.actionable.txt"
     """

     script:

     """
     touch ${rnaseq}.dummy

     Actionable.v3.pl rnaseq ${rnaseq}.dummy ${rnaseq} ${annotation_coding_rare} ${combined_gene_list} ${somatic_actionable_sites} > ${meta.id}.rnaseq.actionable.txt.gl
     Actionable.v3.pl somatic ${somatic_actionable_sites} ${combined_gene_list} ${rnaseq} ${annotation_coding_rare} > ${meta.id}.rnaseq.actionable.txt.som
     germlineOnly.pl ${meta.id}.rnaseq.actionable.txt.gl ${meta.id}.rnaseq.actionable.txt.som > ${meta.id}.rnaseq.actionable.txt
     rm -rf ${meta.id}.rnaseq.actionable.txt.gl ${meta.id}.rnaseq.actionable.txt.som ${rnaseq}.dummy

     """
}


process Actionable_fusion {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/Actionable", mode: "${params.publishDirMode}"

     input:
     path(libs)
     val(meta)
    
     output:
     tuple val(meta),
        path("${meta.id}.fusion.actionable.txt")

     stub:
     """
     touch "${meta.id}.fusion.actionable.txt"
     """

     script:

     """
     if [ \$(wc -l < ${libs.join(' ')}) -ge 2];
     then
        cat ${libs.join(' ')} |sort |uniq > ${meta.id}.fusion.actionable.txt.tmp
        grep "LeftGene" ${meta.id}.fusion.actionable.txt.tmp >${meta.id}.fusion.actionable.txt
        grep -v "LeftGene" ${meta.id}.fusion.actionable.txt.tmp >>${meta.id}.fusion.actionable.txt
        rm -rf ${meta.id}.fusion.actionable.txt.tmp
     else
        touch ${meta.id}.fusion.actionable.txt
     fi

     """
}