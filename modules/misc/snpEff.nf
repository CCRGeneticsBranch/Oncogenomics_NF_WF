process SnpEff {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
        path(vcf),
        path(dbNSFP2_4),
        path(dbNSFP2_4_tbi)

     output:
     tuple val(meta),
        path("${meta.lib}.HC_RNASeq.raw.snpEff.vcf")

     stub:
     """
     touch "${meta.lib}.HC_RNASeq.raw.snpEff.vcf"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """

     set -exo pipefail

     java -jar \$SNPEFF_HOME/SnpSift.jar dbnsfp -db ${dbNSFP2_4}   -a ${vcf} | java -jar \$SNPEFF_HOME/snpEff.jar -t -canon GRCh37.75 > ${prefix}.HC_RNASeq.raw.snpEff.vcf

     """
}


process Vcf2txt {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
        path(vcf)

     output:
     tuple val(meta),
        path("${meta.lib}.HC_RNASeq.snpEff.txt")

     stub:
     """
     touch "${meta.lib}.HC_RNASeq.snpEff.txt"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """

     vcf2txt.pl ${vcf} ./ > ${prefix}.HC_RNASeq.snpEff.txt

     """
}


