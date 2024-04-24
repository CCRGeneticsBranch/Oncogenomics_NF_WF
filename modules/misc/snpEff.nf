process SnpEff {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

     input:
     tuple val(meta),
        path(vcf),
        path(dbNSFP2_4),
        path(dbNSFP2_4_tbi),
        path(Biowulf_snpEff_config),
        val(tool_ch)

     output:
     tuple val(meta),path("${meta.lib}.*${meta.type}.raw.snpEff.vcf"), emit: raw_snpeff
     path "versions.yml"             , emit: versions

     stub:
     """
     touch "${meta.lib}.*${meta.type}.raw.snpEff.vcf"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """

     set -exo pipefail

     java -jar \$SNPEFF_HOME/SnpSift.jar dbnsfp -db ${dbNSFP2_4}  -c ${Biowulf_snpEff_config} -a ${vcf} | java -jar \$SNPEFF_HOME/snpEff.jar -t -canon GRCh37.75 > ${prefix}.${tool_ch}_${meta.type}.raw.snpEff.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Snpsift: \$(java -jar \$SNPEFF_HOME/SnpSift.jar dbnsfp 2>&1 |grep -E '^SnpSift'|cut -f3 -d " ")
    END_VERSIONS
     """
}



process Vcf2txt {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(vcf),val(tool_ch)

     output:
     tuple val(meta),
        path("${meta.lib}.${tool_ch}_${meta.type}.snpEff.txt")

     stub:
     """
     touch "${meta.lib}.${tool_ch}_${meta.type}.snpEff.txt"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """

     vcf2txt.pl ${vcf} ./ > ${prefix}.${tool_ch}_${meta.type}.snpEff.txt
     """
}
