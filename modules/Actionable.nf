process Actionable_variants {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/Actionable", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
      path(dbinput),
      path(annotation_coding_rare),
      path(combined_gene_list),
      path(somatic_actionable_sites),
      val(group)

     output:
     tuple val(meta),
        path("${meta.id}.${group}.actionable.txt")

     stub:
     """
     touch "${meta.id}.${group}.actionable.txt"
     """

     script:

     """
     touch ${dbinput}.dummy

     Actionable.v3.pl rnaseq ${dbinput}.dummy ${dbinput} ${annotation_coding_rare} ${combined_gene_list} ${somatic_actionable_sites} > ${meta.id}.${group}.actionable.txt.gl
     Actionable.v3.pl somatic ${somatic_actionable_sites} ${combined_gene_list} ${dbinput} ${annotation_coding_rare} > ${meta.id}.${group}.actionable.txt.som
     germlineOnly.pl ${meta.id}.${group}.actionable.txt.gl ${meta.id}.${group}.actionable.txt.som > ${meta.id}.${group}.actionable.txt
     rm -rf ${meta.id}.${group}.actionable.txt.gl ${meta.id}.${group}.actionable.txt.som ${dbinput}.dummy

     """
}


process Actionable_fusion {

     tag "$meta.id"

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
     if [ \$(wc -l < ${libs.join(' ')}) -ge 2 ];
     then
        cat ${libs.join(' ')} |sort |uniq > ${meta.id}.fusion.actionable.txt.tmp
        grep "LeftGene" ${meta.id}.fusion.actionable.txt.tmp >${meta.id}.fusion.actionable.txt
        grep -v "LeftGene" ${meta.id}.fusion.actionable.txt.tmp >>${meta.id}.fusion.actionable.txt
        rm -rf ${meta.id}.fusion.actionable.txt.tmp
     else
        echo -e "#LeftGene\tRightGene\tChr_Left\tPosition\tChr_Right\tPosition\tSample\tTool\tSpanReadCount" > ${meta.id}.fusion.actionable.txt
     fi

     """
}