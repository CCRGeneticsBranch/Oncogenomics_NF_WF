process FormatInput {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/annotation", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
        path(vcf2txt),
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
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
     set -exo pipefail
     cat  ${hotspotpileup} |sort > ${prefix}.hotspot
     cut -f 1-5 ${vcf2txt} ${prefix}.hotspot |sort |uniq > AnnotationInput
     MakeAnnotationInputs.pl AnnotationInput

     """
}

process Annovar {
 
     tag "$meta.lib"

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

  tag "$meta.lib"

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
     touch "AnnotationInput"
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

     tag "$meta.lib"

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
        path(snpeff_txt),
        path(ACMG),
        path(hg19_BLsites),
        path(hg19_WLsites)

    output:

   tuple val(meta),path("${meta.id}.Annotations.coding.rare.txt") , emit: rare_annotation       
   tuple val(meta),path("${meta.id}.Annotations.final.txt") , emit: final_annotation
   tuple val(meta),path("${meta.id}.HC_RNASeq.annotated.txt") , emit: hc_RNAseq


     stub:
     """
       touch "${meta.id}.Annotations.coding.rare.txt"
       touch "${meta.id}.Annotations.final.txt"
       touch "${meta.id}.HC_RNASeq.annotated.txt"
     """

   script:
   
     """

     echo "${Anno_input_final}
${gene_out}
${clinseq_out}
${clinvar_out}
${cosmic_out}
${hgmd_out}
${docm_out}
${candl_out}
${tcc_out}
${mcg_out}
${pcg_out}
${civic_out}" > list
     
     CombineAnnotations.pl list > AnnotationInput.annotations.final.txt.tmp     
     GeneAnnotation.pl ${ACMG} AnnotationInput.annotations.final.txt.tmp > ${meta.id}.Annotations.final.txt
     ProteinCodingRare.pl ${hg19_BLsites} ${hg19_WLsites} ${meta.id}.Annotations.final.txt 0.05 > ${meta.id}.Annotations.coding.rare.txt
     addAnnotations2vcf.pl ${meta.id}.Annotations.coding.rare.txt ${snpeff_txt} > ${meta.id}.HC_RNASeq.annotated.txt
     """

}

