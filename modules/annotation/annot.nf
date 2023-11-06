process FormatInput {

     tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/annotation", mode: "${params.publishDirMode}"

     input:
     //tuple val(meta),path(vcf2txt),path(hotspotpileup)
     path(snpeff_vcfs)
     tuple val(meta),path(hotspot)

     output:
     tuple val(meta),
        path("AnnotationInput.anno"),
        path("AnnotationInput.sift"),
        path("AnnotationInput")
     stub:
     """
     touch "AnnotationInput.anno"
     touch "AnnotationInput.sift"
     touch "AnnotationInput"
     """

     script:
     """
     set -exo pipefail
     cat  ${hotspot} |sort > ${meta.id}.hotspot
     cut -f 1-5 ${snpeff_vcfs.join(' ')} ${meta.id}.hotspot |sort |uniq > AnnotationInput
     MakeAnnotationInputs.pl AnnotationInput

     """
}


process FormatInput_TN {

     tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/annotation", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
     path(mutect_snpeff),
     path(strelka_indel_snpeff),
     path(strelka_snvs_snpeff),
     path(hc_normal_snpeff),
     path(hc_tumor_snpeff),
     path(hotspotpileup)


     output:
     tuple val(meta),
        path("AnnotationInput.anno"),
        path("AnnotationInput.sift"),
        path("AnnotationInput")
     stub:
     """
     touch "AnnotationInput.anno"
     touch "AnnotationInput.sift"
     touch "AnnotationInput"
     """

     script:
     """
     set -exo pipefail
     cut -f 1-5 ${mutect_snpeff} ${strelka_indel_snpeff} ${strelka_snvs_snpeff} ${hc_normal_snpeff} ${hc_tumor_snpeff} ${hotspotpileup} |sort |uniq > AnnotationInput
     MakeAnnotationInputs.pl AnnotationInput

     """
}

process Twolib_FormatInput {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/annotation", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(vcf2txt1),path(vcf2txt2)
     tuple val(meta2),path(hotspotpileup)

     output:
     tuple val(meta),
        path("AnnotationInput.anno"),
        path("AnnotationInput.sift"),
        path("AnnotationInput")
     stub:
     """
     touch "AnnotationInput.anno"
     touch "AnnotationInput.sift"
     touch "AnnotationInput"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.id}"
     """
     set -exo pipefail
     cat  ${hotspotpileup} |sort > ${prefix}.hotspot
     cut -f 1-5 ${vcf2txt1} ${vcf2txt2} ${prefix}.hotspot |sort |uniq > AnnotationInput
     MakeAnnotationInputs.pl AnnotationInput

     """
}

process Annovar {

     tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/annotation", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
        path(annotation),
        path(sift),
        path(Anno_input),
        path(annovar_data),
        path(clinseq),
        path(pcg)
     output:
      tuple val(meta),
        path("AnnotationInput.cosmic"),
        path("AnnotationInput.clinseq"),
        path("AnnotationInput.cadd"),
        path("AnnotationInput.pcg"),
        path("AnnotationInput.gene")

     stub:
     """
     touch "AnnotationInput.cosmic"
     touch "AnnotationInput.clinseq"
     touch "AnnotationInput.cadd"
     touch "AnnotationInput.pcg"
     touch "AnnotationInput.gene"
     """

     script:

     """
     set -exo pipefail
     table_annovar.pl ${annotation} ${annovar_data} -buildver hg19 --dot2underline -out AnnotationInput.anno.cosmic -remove -protocol cosmic78 -operation f -polish -nastring "NA"
     mv AnnotationInput.anno.cosmic.hg19_multianno.txt AnnotationInput.cosmic

     annotate_variation.pl  ${annotation} ${annovar_data} -buildver hg19 -otherinfo --outfile AnnotationInput.anno.clinseq -filter -dbtype generic -genericdbfile `basename ${clinseq} `
     awk '{OFS="\t"};{print \$3,\$4,\$5,\$6,\$7,\$2}' AnnotationInput.anno.clinseq.hg19_generic_dropped |sed -e 's/,/\t/g' > AnnotationInput.clinseq
     head -1 ${annovar_data}/hg19_clinseq_951.txt >>AnnotationInput.clinseq

     annotate_variation.pl ${annotation} ${annovar_data} -buildver hg19 -otherinfo  -filter -dbtype cadd
     annotate_variation.pl AnnotationInput.anno.hg19_cadd_filtered ${annovar_data} -buildver hg19 -otherinfo -filter -dbtype caddindel
     cut -f 2-7 AnnotationInput.anno.hg19_cadd_dropped AnnotationInput.anno.hg19_cadd_filtered.hg19_caddindel_dropped|sed -e 's/,/\t/g' |awk '{OFS="\t"};{print \$3,\$4,\$5,\$6,\$7,\$1,\$2}' >AnnotationInput.cadd
     head -1 ${annovar_data}/hg19_caddindel.txt >>AnnotationInput.cadd

     annotate_variation.pl ${annotation} ${annovar_data} -buildver hg19 -otherinfo --outfile AnnotationInput.anno.pcg -filter -dbtype generic -genericdbfile `basename ${pcg} `
     awk -F "\t" '{OFS="\t"};{print \$3,\$4,\$5,\$6,\$7,\$2}' AnnotationInput.anno.pcg.hg19_generic_dropped |sed -e 's/,/\t/g' >AnnotationInput.pcg
     head -1 ${annovar_data}/hg19_PCG_042616.txt >>AnnotationInput.pcg

     table_annovar.pl  ${annotation} ${annovar_data} -buildver hg19 -out AnnotationInput.anno.gene -remove -protocol refGene,cytoBand,snp138,1000g2014oct_all,1000g2014oct_eur,1000g2014oct_afr,1000g2014oct_amr,1000g2014oct_eas,1000g2014oct_sas,esp6500_all,esp6500_ea,esp6500_aa,exac03nontcga,exac03,cg69,nci60 -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring "-1" -polish  --argument "-hgvs",,,,,,,,,,,,,,,
     mv AnnotationInput.anno.gene.hg19_multianno.txt AnnotationInput.gene

     sed -i '1s/[.]/_/g' AnnotationInput.gene

     """
}

process Custom_annotation {

  tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/annotation", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
        path(annotation),
        path(sift),
        path(Anno_input),
        path(clinvar),
        path(hgmd),
        path(matchTrial),
        path(mcg),
        path(DoCM),
        path(CanDL),
        path(targetted_cancer_care),
        path(civic)
     output:
     tuple val(meta),
        path("AnnotationInput.clinvar"),
        path("AnnotationInput.hgmd"),
        path("AnnotationInput.match"),
        path("AnnotationInput.mcg"),
        path("AnnotationInput.docm"),
        path("AnnotationInput.candl"),
        path("AnnotationInput.tcc"),
        path("AnnotationInput.civic"),
        path("AnnotationInput_final")

     stub:
     """
     touch "AnnotationInput.clinvar"
     touch "AnnotationInput.hgmd"
     touch "AnnotationInput.match"
     touch "AnnotationInput.mcg"
     touch "AnnotationInput.docm"
     touch "AnnotationInput.candl"
     touch "AnnotationInput.tcc"
     touch "AnnotationInput.civic"
     touch "AnnotationInput_final"
     """

     script:

     """
       addAnnotation.pl ${clinvar} ${annotation} >AnnotationInput.clinvar
       addAnnotation.pl ${hgmd} ${annotation} >AnnotationInput.hgmd
       addAnnotation.pl ${matchTrial} ${annotation} >AnnotationInput.match
       addAnnotation.pl ${mcg} ${annotation} >AnnotationInput.mcg
       addAnnotation.pl ${DoCM} ${annotation} >AnnotationInput.docm
       addAnnotation.pl ${CanDL} ${annotation} >AnnotationInput.candl
       addAnnotation.pl	${targetted_cancer_care} ${annotation} >AnnotationInput.tcc
       addAnnotation.pl	${civic} ${annotation} >AnnotationInput.civic
       cp ${Anno_input} AnnotationInput_final
     """
}


process Combine_annotation {

     tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/annotation", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
        path(cosmic_out),
        path(clinseq_out),
        path(cadd_out),
        path(pcg_out),
        path(gene_out),
       	path(clinvar_out),
       	path(hgmd_out),
       	path(match_out),
       	path(mcg_out),
       	path(docm_out),
       	path(candl_out),
       	path(tcc_out),
       	path(civic_out),
        path(Anno_input_final),
        path(ACMG),
        path(hg19_BLsites),
        path(hg19_WLsites)

    output:

   tuple val(meta),path("${meta.id}.Annotations.coding.rare.txt") , emit: rare_annotation
   tuple val(meta),path("${meta.id}.Annotations.final.txt") , emit: final_annotation

     stub:
     """
       touch "${meta.id}.Annotations.coding.rare.txt"
       touch "${meta.id}.Annotations.final.txt"
     """

   script:

     """

echo "Anno_input_final: ${Anno_input_final}"
echo "gene_out: ${gene_out}"
echo "clinseq_out: ${clinseq_out}"
echo "cadd_out: ${cadd_out}"
echo "clinvar_out: ${clinvar_out}"
echo "cosmic_out: ${cosmic_out}"
echo "hgmd_out: ${hgmd_out}"
echo "match_out: ${match_out}"
echo "docm_out: ${docm_out}"
echo "candl_out: ${candl_out}"
echo "tcc_out: ${tcc_out}"
echo "mcg_out: ${mcg_out}"
echo "pcg_out: ${pcg_out}"
echo "civic_out: ${civic_out}"

echo "${Anno_input_final}
${gene_out}
${clinseq_out}
${cadd_out}
${clinvar_out}
${cosmic_out}
${hgmd_out}
${match_out}
${docm_out}
${candl_out}
${tcc_out}
${mcg_out}
${civic_out}
${pcg_out}" > list

     CombineAnnotations.pl list > AnnotationInput.annotations.final.txt.tmp
     GeneAnnotation.pl ${ACMG} AnnotationInput.annotations.final.txt.tmp > ${meta.id}.Annotations.final.txt
     ProteinCodingRare.pl ${hg19_BLsites} ${hg19_WLsites} ${meta.id}.Annotations.final.txt 0.05 > ${meta.id}.Annotations.coding.rare.txt
     """

}

process AddAnnotation {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(snpeff_txt),path(rare_annotation)


     output:
     tuple val(meta),path("${meta.lib}.HC_${meta.type}.annotated.txt") , emit: hc_anno_txt

     stub:
     """
       touch "${meta.lib}.HC_${meta.type}.annotated.txt"
     """
     script:

     """

     addAnnotations2vcf.pl  ${rare_annotation} ${snpeff_txt}  > ${meta.lib}.HC_${meta.type}.annotated.txt

     """

}


process AddAnnotation_TN {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(N_snpeff_txt),path(T_snpeff_txt),path(rare_annotation)


     output:
     tuple val(meta),path("${meta.lib}.HC_${meta.type}.annotated.txt") , emit: Tumor_hc_anno_txt
     tuple val(meta),path("${meta.normal_id}.HC_${meta.normal_type}.annotated.txt") , emit: Normal_hc_anno_txt
     stub:
     """
       touch "${meta.lib}.HC_${meta.type}.annotated.txt"
       touch "${meta.normal_id}.HC_${meta.normal_type}.annotated.txt"
     """
     script:

     """

     addAnnotations2vcf.pl  ${rare_annotation} ${T_snpeff_txt}  > ${meta.lib}.HC_${meta.type}.annotated.txt
     addAnnotations2vcf.pl  ${rare_annotation} ${N_snpeff_txt}  > ${meta.normal_id}.HC_${meta.normal_type}.annotated.txt

     """

}

process AddAnnotation_somatic_variants {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
     path(mutect_txt),
     path(strelka_indels_txt),
     path(strelka_snvs_txt),
     path(rare_annotation)

     output:
     tuple val(meta),
     path("${meta.lib}.MuTect.annotated.txt"),
     path("${meta.lib}.strelka.indels.annotated.txt"),
     path("${meta.lib}.strelka.snvs.annotated.txt")

     stub:
     """
       touch "${meta.lib}.MuTect.annotated.txt"
       touch "${meta.lib}.strelka.indels.annotated.txt"
       touch "${meta.lib}.strelka.snvs.annotated.txt"

     """
     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """

     addAnnotations2vcf.pl  ${rare_annotation} ${mutect_txt}  > ${prefix}.MuTect.annotated.txt
     addAnnotations2vcf.pl  ${rare_annotation} ${strelka_indels_txt}  > ${prefix}.strelka.indels.annotated.txt
     addAnnotations2vcf.pl  ${rare_annotation} ${strelka_snvs_txt}  > ${prefix}.strelka.snvs.annotated.txt

    """

}


process AddAnnotationFull_somatic_variants {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
     path(mutect_txt),
     path(strelka_indels_txt),
     path(strelka_snvs_txt),
     path(final_annotation)

     output:
     tuple val(meta),
     path("${meta.lib}.MuTect.annotatedFull.txt"),
     path("${meta.lib}.strelka.indels.annotatedFull.txt"),
     path("${meta.lib}.strelka.snvs.annotatedFull.txt")


     stub:
     """
       touch "${meta.lib}.MuTect.annotatedFull.txt"
       touch "${meta.lib}.strelka.indels.annotatedFull.txt"
       touch "${meta.lib}.strelka.snvs.annotatedFull.txt"

     """
     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
     addAnnotations2vcf.pl  ${final_annotation} ${mutect_txt}  > ${prefix}.MuTect.annotatedFull.txt
     addAnnotations2vcf.pl  ${final_annotation} ${strelka_indels_txt}  > ${prefix}.strelka.indels.annotatedFull.txt
     addAnnotations2vcf.pl  ${final_annotation} ${strelka_snvs_txt}  > ${prefix}.strelka.snvs.annotatedFull.txt

     """

}
