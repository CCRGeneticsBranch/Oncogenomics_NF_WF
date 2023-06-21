process SnpEff {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
        path(vcf),
        path(dbNSFP2_4),
        path(dbNSFP2_4_tbi),
        path(Biowulf_snpEff_config)

     output:
     tuple val(meta),
        path("${meta.lib}.HC_${meta.type}.raw.snpEff.vcf")

     stub:
     """
     touch "${meta.lib}.HC_${meta.type}.raw.snpEff.vcf"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """

     set -exo pipefail

     java -jar \$SNPEFF_HOME/SnpSift.jar dbnsfp -db ${dbNSFP2_4}  -c ${Biowulf_snpEff_config} -a ${vcf} | java -jar \$SNPEFF_HOME/snpEff.jar -t -canon GRCh37.75 > ${prefix}.HC_${meta.type}.raw.snpEff.vcf

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
        path("${meta.lib}.HC_${meta.type}.snpEff.txt")

     stub:
     """
     touch "${meta.lib}.HC_${meta.type}.snpEff.txt"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """

     vcf2txt.pl ${vcf} ./ > ${prefix}.HC_${meta.type}.snpEff.txt

     """
}


