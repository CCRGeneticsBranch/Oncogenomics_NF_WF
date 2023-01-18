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
        path("AnnotationInput.sift")

     stub:
     """
     touch "AnnotationInput.anno"
     touch "AnnotationInput.sift"
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
     sed -i '1s|.|_|g' AnnotationInput.gene
     '''
}
