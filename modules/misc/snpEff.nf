process SnpEff {

     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/GATK_RNAseq", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
        path(vcf),
        path(dbNSFP2_4),
        path(dbNSFP2_4_tbi)

     output:
     tuple val("${dataset_id}"),
        path("${dataset_id}.HC_RNASeq.raw.snpEff.vcf")

     stub:
     """
     touch "${dataset_id}.HC_RNASeq.raw.snpEff.vcf"
     """

     shell:
     '''
     set -exo pipefail

     java -jar \$SNPEFF_HOME/SnpSift.jar dbnsfp -db !{dbNSFP2_4}   -a !{vcf} | java -jar \$SNPEFF_HOME/snpEff.jar -t -canon GRCh37.75 > !{dataset_id}.HC_RNASeq.raw.snpEff.vcf

     '''
}


process Vcf2txt {

     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/GATK_RNAseq", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
        path(vcf)

     output:
     tuple val("${dataset_id}"),
        path("${dataset_id}.HC_RNASeq.snpEff.txt")

     stub:
     """
     touch "${dataset_id}.HC_RNASeq.snpEff.txt"
     """

     shell:
     '''
     set -exo pipefail

     vcf2txt.pl !{vcf} ./ > !{dataset_id}.HC_RNASeq.snpEff.txt

     '''
}


