process SnpEff {

     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/${library}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
        val(library),
        path(vcf),
        path(dbNSFP2_4),
        path(dbNSFP2_4_tbi)

     output:
     tuple val("${dataset_id}"),
        val("$library"),
        path("${library}.HC_RNASeq.raw.snpEff.vcf")

     stub:
     """
     touch "${library}.HC_RNASeq.raw.snpEff.vcf"
     """

     shell:
     '''
     set -exo pipefail

     java -jar \$SNPEFF_HOME/SnpSift.jar dbnsfp -db !{dbNSFP2_4}   -a !{vcf} | java -jar \$SNPEFF_HOME/snpEff.jar -t -canon GRCh37.75 > !{library}.HC_RNASeq.raw.snpEff.vcf

     '''
}


process Vcf2txt {

     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/${library}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
        val(library),
        path(vcf)

     output:
     tuple val("${dataset_id}"),
        val("$library"),
        path("${library}.HC_RNASeq.snpEff.txt")

     stub:
     """
     touch "${library}.HC_RNASeq.snpEff.txt"
     """

     shell:
     '''
     set -exo pipefail

     vcf2txt.pl !{vcf} ./ > !{library}.HC_RNASeq.snpEff.txt

     '''
}


