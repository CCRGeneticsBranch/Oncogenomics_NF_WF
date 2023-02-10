process FormatInput {

     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/annotation", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
        path(vcf2txt),
        path(hotspotpileup)

     output:
     tuple val(dataset_id),
        path("AnnotationInput.anno"),
        path("AnnotationInput.sift"),
        path("AnnotationInput")
     stub:
     """
     touch "AnnotationInput.anno"
     touch "AnnotationInput.sift"
     touch "AnnotationInput"
     """

     shell:
     '''
     set -exo pipefail
     cat  !{hotspotpileup} |sort > !{dataset_id}.hotspot
     cut -f 1-5 !{vcf2txt} !{dataset_id}.hotspot |sort |uniq > AnnotationInput
     MakeAnnotationInputs.pl AnnotationInput

     '''
}

process Annovar {
 
     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/annotation", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
        path(annotation),
        path(sift),
        path(Anno_input),
        path(annovar_data),
        path(clinseq),
        path(pcg)
     output:
     tuple val(dataset_id),
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
        
     shell:
     '''
     set -exo pipefail
     table_annovar.pl !{annotation} !{annovar_data} -buildver hg19 --dot2underline -out AnnotationInput.anno.cosmic -remove -protocol cosmic78 -operation f -polish -nastring "NA"
     mv AnnotationInput.anno.cosmic.hg19_multianno.txt AnnotationInput.cosmic
     annotate_variation.pl  !{annotation} !{annovar_data} -buildver hg19 -otherinfo --outfile AnnotationInput.anno.clinseq -filter -dbtype generic -genericdbfile `basename !{clinseq} `
     awk '{OFS="\t"};{print $3,$4,$5,$6,$7,$2}' AnnotationInput.anno.clinseq.hg19_generic_dropped |sed -e 's/,/\t/g' > AnnotationInput.clinseq
     head -1 !{annovar_data}/hg19_clinseq_951.txt >>AnnotationInput.clinseq
     annotate_variation.pl !{annotation} !{annovar_data} -buildver hg19 -otherinfo  -filter -dbtype cadd
     annotate_variation.pl AnnotationInput.anno.hg19_cadd_filtered !{annovar_data} -buildver hg19 -otherinfo -filter -dbtype caddindel
     cut -f 2-7 AnnotationInput.anno.hg19_cadd_dropped AnnotationInput.anno.hg19_cadd_filtered.hg19_caddindel_dropped|sed -e 's/,/\t/g' |awk '{OFS="\t"};{print $3,$4,$5,$6,$7,$1,$2}' >AnnotationInput.cadd
     head -1 !{annovar_data}/hg19_caddindel.txt >>AnnotationInput.cadd
     annotate_variation.pl !{annotation} !{annovar_data} -buildver hg19 -otherinfo --outfile AnnotationInput.anno.pcg -filter -dbtype generic -genericdbfile `basename !{pcg} `
     awk -F "\t" '{OFS="\t"};{print $3,$4,$5,$6,$7,$2}' AnnotationInput.anno.pcg.hg19_generic_dropped |sed -e 's/,/\t/g' >AnnotationInput.pcg
     head -1 !{annovar_data}/hg19_PCG_042616.txt >>AnnotationInput.pcg
     table_annovar.pl  !{annotation} !{annovar_data} -buildver hg19 -out AnnotationInput.anno.gene -remove -protocol refGene,cytoBand,snp138,1000g2014oct_all,1000g2014oct_eur,1000g2014oct_afr,1000g2014oct_amr,1000g2014oct_eas,1000g2014oct_sas,esp6500_all,esp6500_ea,esp6500_aa,exac03nontcga,exac03,cg69,nci60 -operation g,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f -nastring "-1" -polish  --argument "-hgvs",,,,,,,,,,,,,,,
     mv AnnotationInput.anno.gene.hg19_multianno.txt AnnotationInput.gene
#     sed -i '1s|.|_|g' AnnotationInput.gene
     sed -i '1s/[.]/_/g' AnnotationInput.gene

     '''
}

process Custom_annotation {

     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/annotation", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
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
     tuple val(dataset_id),
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

     shell:
     '''
       addAnnotation.pl !{clinvar} !{annotation} >AnnotationInput.clinvar
       addAnnotation.pl !{hgmd} !{annotation} >AnnotationInput.hgmd
       addAnnotation.pl !{matchTrial} !{annotation} >AnnotationInput.match
       addAnnotation.pl !{mcg} !{annotation} >AnnotationInput.mcg
       addAnnotation.pl !{DoCM} !{annotation} >AnnotationInput.docm
       addAnnotation.pl !{CanDL} !{annotation} >AnnotationInput.candl
       addAnnotation.pl	!{targetted_cancer_care} !{annotation} >AnnotationInput.tcc
       addAnnotation.pl	!{civic} !{annotation} >AnnotationInput.civic
       cp !{Anno_input} AnnotationInput_final
     '''
}


process Combine_annotation {

     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/annotation", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
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
    tuple val("${dataset_id}"),
        path("${dataset_id}.Annotations.coding.rare.txt")
        path("${dataset_id}.Annotations.final.txt")
        path("${dataset_id}.HC_RNASeq.annotated.txt")
        
     stub:
     """
       touch "${dataset_id}.Annotations.coding.rare.txt"
       touch "${dataset_id}.Annotations.final.txt"
       touch "${dataset_id}.HC_RNASeq.annotated.txt"
     """

     shell:
     '''
     echo "!{Anno_input_final}
!{gene_out}
!{clinseq_out}
!{cadd_out}
!{clinvar_out}
!{cosmic_out}
!{hgmd_out}
!{match_out}
!{docm_out}
!{candl_out}
!{tcc_out}
!{mcg_out}
!{pcg_out}
!{civic_out}" > list
     
     CombineAnnotations.pl list > AnnotationInput.annotations.final.txt.tmp     
     GeneAnnotation.pl !{ACMG} AnnotationInput.annotations.final.txt.tmp > !{dataset_id}.Annotations.final.txt
     ProteinCodingRare.pl !{hg19_BLsites} !{hg19_WLsites} !{dataset_id}.Annotations.final.txt 0.05 > !{dataset_id}.Annotations.coding.rare.txt
     addAnnotations2vcf.pl !{dataset_id}.Annotations.coding.rare.txt !{snpeff_txt} > !{dataset_id}.HC_RNASeq.annotated.txt
     '''
}

