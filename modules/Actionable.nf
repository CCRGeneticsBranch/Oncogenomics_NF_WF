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