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
     if [[ "${meta.type}" == "tumor_RNA" || "${meta.type}" == "cell_line_RNA" || "${meta.type}" == "xeno_RNA" ]]; then
          tag="rnaseq"
          makeDBVariantFile.pl ${dbinput_annot_libs.join(' ')}|sed 's/trim_//'|AddSampleType.pl - "${meta.lib} ${meta.type}" "${meta.lib} ${meta.sc}" > ${meta.id}.\$tag

     elif [[ "${meta.type}" == "normal_DNA" || "${meta.type}" == "blood_DNA" ]]; then
          tag="germline"
          makeDBVariantFile.pl ${dbinput_annot_libs.join(' ')}|sed 's/trim_//'|AddSampleType.pl - "${meta.lib} ${meta.type}" "${meta.lib} ${meta.sc}" > ${meta.id}.\$tag.tmp
          addFS.pl ${meta.id}.\$tag.tmp ${dbinput_snpeff_libs.join(' ')} >${meta.id}.\$tag
          rm -rf ${meta.id}.\$tag.tmp
     elif [[ ( "${meta.type}" == "tumor_DNA" || "${meta.type}" == "cell_line_DNA" ) && "${meta.normal_type}" != "normal_DNA" && "${meta.normal_type}" != "blood_DNA" && "${meta.rna_type}" != "tumor_RNA" && "${meta.rna_type}" != "cell_line_RNA" && "${meta.rna_type}" != "xeno_RNA" ]]; then
          tag="variants"
          makeDBVariantFile.pl ${dbinput_annot_libs.join(' ')}|sed 's/trim_//'|AddSampleType.pl - "${meta.lib} ${meta.type}" "${meta.lib} ${meta.sc}" > ${meta.id}.\$tag.tmp
          addFS.pl ${meta.id}.\$tag.tmp ${dbinput_snpeff_libs.join(' ')} >${meta.id}.\$tag
          rm -rf ${meta.id}.\$tag.tmp
     elif [[ "${meta.type}" == "tumor_DNA" && ( "${meta.normal_type}" == "normal_DNA" || "${meta.normal_type}" == "blood_DNA" ) && "${meta.rna_type}" != "tumor_RNA" && "${meta.rna_type}" != "xeno_RNA" ]]; then
          tag1="somatic"
          tag2="germline"
     elif [[ ( "${meta.type}" == "tumor_DNA" || "${meta.type}" == "cell_line_DNA" ) && ( "${meta.rna_type}" == "tumor_RNA" || "${meta.rna_type}" == "cell_line_RNA" ) && "${meta.normal_type}" != "normal_DNA" && "${meta.normal_type}" != "blood_DNA" ]]; then
          tag1="variants"
          tag2="rnaseq"
          makeDBVariantFile.pl ${dbinput_annot_libs[0]}|sed 's/trim_//'|AddSampleType.pl - "${meta.rna_lib} ${meta.rna_type} ${meta.lib} ${meta.type}" "${meta.rna_lib} ${meta.RNA_sc} ${meta.lib} ${meta.sc}" > ${meta.id}.\$tag1.tmp
          addFS.pl ${meta.id}.\$tag1.tmp ${dbinput_snpeff_libs.join(' ')} >${meta.id}.\$tag1
          makeDBVariantFile.pl ${dbinput_annot_libs[1]} |sed 's/trim_//'|AddSampleType.pl - "${meta.rna_lib} ${meta.rna_type} ${meta.lib} ${meta.type}" "${meta.rna_lib} ${meta.RNA_sc} ${meta.lib} ${meta.sc}" > ${meta.id}.\$tag2
          rm *.tmp
     elif [[ ( "${meta.normal_type}" == "normal_DNA" || "${meta.normal_type}" == "blood_DNA" ) && "${meta.type}" == "tumor_DNA" && ( "${meta.rna_type}" == "tumor_RNA" || "${meta.rna_type}" == "xeno_RNA" ) ]]; then
          tag1="somatic"
          tag2="germline"
          tag3="rnaseq"
          #somatic
          makeDBVariantFile.pl ${dbinput_annot_libs[3]} ${dbinput_annot_libs[5]} ${dbinput_annot_libs[4]} |sed 's/trim_//'|AddSampleType.pl - "${meta.normal_id} ${meta.normal_type} ${meta.lib} ${meta.type} ${meta.rna_lib} ${meta.rna_type}" "${meta.normal_id} ${meta.N_sc} ${meta.lib} ${meta.sc} ${meta.rna_lib} ${meta.RNA_sc}" > ${meta.id}.\$tag1
          #germline
          makeDBVariantFile.pl ${dbinput_annot_libs[0]} ${dbinput_annot_libs[1]} |sed 's/trim_//'|AddSampleType.pl - "${meta.normal_id} ${meta.normal_type} ${meta.lib} ${meta.type} ${meta.rna_lib} ${meta.rna_type}" "${meta.normal_id} ${meta.N_sc} ${meta.lib} ${meta.sc} ${meta.rna_lib} ${meta.RNA_sc}" > ${meta.id}.\$tag2.tmp
          addFS.pl ${meta.id}.\$tag2.tmp ${dbinput_snpeff_libs.join(' ')} >${meta.id}.\$tag2
          #rnaseq
          makeDBVariantFile.pl ${dbinput_annot_libs[2]} |sed 's/trim_//'|AddSampleType.pl - "${meta.normal_id} ${meta.normal_type} ${meta.lib} ${meta.type} ${meta.rna_lib} ${meta.rna_type}" "${meta.normal_id} ${meta.N_sc} ${meta.lib} ${meta.sc} ${meta.rna_lib} ${meta.RNA_sc}" > ${meta.id}.\$tag3
          rm *.tmp

     fi

     """
}



process DBinput_multiple {

     tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.id}/db", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(annotated), path(snpeff)

     output:
     tuple val(meta),
        path("${meta.id}*")

     stub:
     """
     touch "${meta.id}*"
     """

     script:

     """
     LIB1=`echo ${annotated[0]}|awk -F"." '{print \$1}'`
     LIB2=`echo ${annotated[1]}|awk -F"." '{print \$1}'`

     if [[ "${meta.type}" == *"RNA"* ]]; then
          tag="rnaseq"
          makeDBVariantFile.pl ${annotated.join(' ')}|sed 's/trim_//'|AddSampleType.pl - "\$LIB1 ${meta.type} \$LIB2 ${meta.type}" "\$LIB1 ${meta.sc} \$LIB2 ${meta.sc}" > ${meta.id}.\$tag.tmp
          mv ${meta.id}.\$tag.tmp ${meta.id}.\$tag
     else
          tag="variants"
          makeDBVariantFile.pl ${annotated.join(' ')}|sed 's/trim_//'|AddSampleType.pl - "\$LIB1 ${meta.type} \$LIB2 ${meta.type}" "\$LIB1 ${meta.sc} \$LIB2 ${meta.sc}" > ${meta.id}.\$tag.tmp
          addFS.pl ${meta.id}.\$tag.tmp ${snpeff.join(' ')} > ${meta.id}.\$tag
     fi

     """
}

process DBinput_multiple_new {

     tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.id}/db", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),val(libtype), val(libsc), path(annot)
     path(snpeff)

     output:
     tuple val(meta),
        path("${meta.id}*")

     stub:
     """
     touch "${meta.id}*"
     """

     script:

     """

     if [[ "${libtype}" == *"RNA"* ]]; then
          tag="rnaseq"
          makeDBVariantFile.pl ${annot.join(' ')}|sed 's/trim_//'|AddSampleType.pl - "${libtype.join(' ')}" "${libsc.join(' ')}" > ${meta.id}.\$tag.tmp
          mv ${meta.id}.\$tag.tmp ${meta.id}.\$tag
     else
          tag="variants"
          makeDBVariantFile.pl ${annot.join(' ')}|sed 's/trim_//'|AddSampleType.pl - "${libtype.join(' ')}" "${libsc.join(' ')}" > ${meta.id}.\$tag.tmp
          addFS.pl ${meta.id}.\$tag.tmp ${snpeff.join(' ')} > ${meta.id}.\$tag
     fi

     """
}


process DBinput_exome_rnaseq {

     tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.id}/db", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),val(libtype), val(libsc), path(dna), path (rna)
     path(snpeff)

     output:
     tuple val(meta),
        path("${meta.id}*")

     stub:
     """
     touch "${meta.id}*"
     """

     script:

     """
     makeDBVariantFile.pl ${rna.join(' ')}|sed 's/trim_//'|AddSampleType.pl - "${libtype.join(' ')}" "${libsc.join(' ')}" > ${meta.id}.rnaseq

     makeDBVariantFile.pl ${dna.join(' ')}|sed 's/trim_//'|AddSampleType.pl - "${libtype.join(' ')}" "${libsc.join(' ')}" > ${meta.id}.variants.tmp
     addFS.pl ${meta.id}.variants.tmp ${snpeff.join(' ')} > ${meta.id}.variants
     rm ${meta.id}.variants.tmp

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
     makeDBVariantFile.pl ${annot_mutect} ${annot_strelka1} ${annot_strelka2}|sed 's/trim_//'|AddSampleType.pl - "${meta.normal_id} ${meta.normal_type} ${meta.lib} ${meta.type}" "${meta.normal_id} ${meta.N_sc} ${meta.lib} ${meta.sc}" > ${meta.id}.${somatic}

     ##Germline
     makeDBVariantFile.pl ${annot_tumor} ${annot_normal} | AddSampleType.pl - "${meta.normal_id} ${meta.normal_type} ${meta.lib} ${meta.type}" "${meta.normal_id} ${meta.N_sc} ${meta.lib} ${meta.sc}" > ${meta.id}.${germline}.tmp
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
