process DBinput {

     tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.id}/db", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(dbinput_annot_libs),path(dbinput_snpeff_libs)

     output:
     tuple val(meta),path("${meta.id}.*")

     stub:
     """
     touch "${meta.id}.*"
     """

     script:

     """
     if [[ "${meta.type}" == "tumor_RNA" || "${meta.type}" == "cell_line_RNA" ]]; then
          tag="rnaseq"
          makeDBVariantFile.pl ${dbinput_annot_libs.join(' ')}|sed 's/trim_//'|AddSampleType.pl - "${meta.lib} ${meta.type}" "${meta.lib} ${meta.sc}" > ${meta.id}.\$tag.tmp
          addFS.pl ${meta.id}.\$tag.tmp ${dbinput_snpeff_libs.join(' ')} >${meta.id}.\$tag
          rm -rf ${meta.id}.\$tag.tmp

     elif [[ "${meta.type}" == "normal_DNA" ]]; then
          tag="germline"
          makeDBVariantFile.pl ${dbinput_annot_libs.join(' ')}|sed 's/trim_//'|AddSampleType.pl - "${meta.lib} ${meta.type}" "${meta.lib} ${meta.sc}" > ${meta.id}.\$tag.tmp
          addFS.pl ${meta.id}.\$tag.tmp ${dbinput_snpeff_libs.join(' ')} >${meta.id}.\$tag
          rm -rf ${meta.id}.\$tag.tmp

     elif [[ "${meta.type}" == "tumor_DNA" || "${meta.type}" == "cell_line_DNA" && "${meta.normal_type}" == "null" ]]; then
          tag="variants"
          makeDBVariantFile.pl ${dbinput_annot_libs.join(' ')}|sed 's/trim_//'|AddSampleType.pl - "${meta.lib} ${meta.type}" "${meta.lib} ${meta.sc}" > ${meta.id}.\$tag.tmp
          addFS.pl ${meta.id}.\$tag.tmp ${dbinput_snpeff_libs.join(' ')} >${meta.id}.\$tag
          rm -rf ${meta.id}.\$tag.tmp

     elif [[ "${meta.type}" == "tumor_DNA" && "${meta.normal_type}" == "normal_DNA" && "${meta.rna_type}" == "null" ]]; then
          tag1="somatic"
          tag2="germline"
     elif [[ "${meta.type}" == "tumor_DNA" && "${meta.normal_type}" == "normal_DNA" && "${meta.rna_type}" == "tumor_RNA" ]]; then
          tag1="somatic"
          tag2="germline"
          tag3="rnaseq"
     fi

     """
}



process DBinput_multiple {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.id}/db", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(annotated1), path(annotated2)
     tuple val(meta),path(snpeff1), path(snpeff2)



     output:
     tuple val(meta),
        path("${meta.id}.${meta.type}")

     stub:
     """
     touch "${meta.id}.${meta.type}"
     """

     script:

     """
     LIB1=`echo ${annotated1}|awk -F"." '{print \$1}'`
     LIB2=`echo ${annotated2}|awk -F"." '{print \$1}'`
    makeDBVariantFile.pl ${annotated1} ${annotated2}|sed 's/trim_//'|AddSampleType.pl - "\$LIB1 ${meta.type} \$LIB2 ${meta.type}" "\$LIB1 ${meta.sc} \$LIB2 ${meta.sc}" > ${meta.id}.${meta.type}.tmp
    if [ ${meta.type}  == 'germline' ] || [ ${meta.type} == 'variants' ]; then
         addFS.pl ${meta.id}.${meta.type}.tmp ${snpeff1} ${snpeff2} >${meta.id}.${meta.type}
         rm -rf ${meta.id}.${meta.type}.tmp
    else
        mv ${meta.id}.${meta.type}.tmp ${meta.id}.${meta.type}
    fi
     """
}

process DBinput_multiples {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.id}/db", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
     path(annot_mutect),path(annot_strelka1),
     path(annot_strelka2),path(somatic_mutect_snpeff),path(somatic_strelka1_snpeff),
     path(somatic_strelka2_snpeff),path(HC_snpeff1),path(HC_snpeff2),
     path(annot_tumor),path(annot_normal),
     val(somatic),val(germline)

     output:
     tuple val(meta),path("${meta.id}.${somatic}")
     tuple val(meta),path("${meta.id}.${germline}")

     stub:
     """
     touch "${meta.id}.${somatic}"
     touch "${meta.id}.${germline}"
     """

     script:

     """
     ##Somatic
     makeDBVariantFile.pl ${annot_mutect} ${annot_strelka1} ${annot_strelka2}|sed 's/trim_//'|AddSampleType.pl - "${meta.normal_id} ${meta.normal_type} ${meta.lib} ${meta.type}" "${meta.normal_id} ${meta.N_sc} ${meta.lib} ${meta.T_sc}" > ${meta.id}.${somatic}.tmp
     addFS.pl ${meta.id}.${somatic}.tmp ${HC_snpeff1} ${HC_snpeff2} ${somatic_mutect_snpeff} ${somatic_strelka1_snpeff} ${somatic_strelka2_snpeff} > ${meta.id}.${somatic}
     rm -rf ${meta.id}.${somatic}.tmp

     ##Germline
     makeDBVariantFile.pl ${annot_tumor} ${annot_normal} | AddSampleType.pl - "${meta.normal_id} ${meta.normal_type} ${meta.lib} ${meta.type}" "${meta.normal_id} ${meta.N_sc} ${meta.lib} ${meta.T_sc}" > ${meta.id}.${germline}.tmp
     addFS.pl ${meta.id}.${germline}.tmp ${HC_snpeff1} ${HC_snpeff2} ${somatic_mutect_snpeff} ${somatic_strelka1_snpeff} ${somatic_strelka2_snpeff} > ${meta.id}.${germline}
     rm -rf ${meta.id}.${germline}.tmp

     """
}


process Somatic_actionable {

     tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/Actionable", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(annotation_coding_rare)
     tuple val(meta2),path(somatic)
     path somatic_actionable_sites
     path combined_gene_list

     output:
     tuple val(meta),
        path("${meta.id}.somatic.actionable.txt")

     stub:
     """
     touch "${meta.id}.somatic.actionable.txt"
     """

     script:
     """
     Actionable.v3.pl somatic ${somatic_actionable_sites} ${combined_gene_list} ${somatic} ${annotation_coding_rare} > ${meta.id}.somatic.actionable.txt

     """
}

process Germline_actionable {

     tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/Actionable", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(annotation_coding_rare)
     tuple val(meta2),path(germline)
     path somatic_actionable_sites
     path combined_gene_list
     tuple val(meta2),path(dbinput_somatic)

     output:
     tuple val(meta),
        path("${meta.id}.germline.actionable.txt")

     stub:
     """
     touch "${meta.id}.germline.actionable.txt"
     """

     script:
     """
     Actionable.v3.pl germline ${dbinput_somatic} ${germline} ${annotation_coding_rare} ${combined_gene_list} ${somatic_actionable_sites}  > ${meta.id}.germline.actionable.txt

     """
}
