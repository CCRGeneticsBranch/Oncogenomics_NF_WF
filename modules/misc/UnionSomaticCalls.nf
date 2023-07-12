process UnionSomaticCalls {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/Actionable", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(mutect_fulltxt), path(strelka_indels_fulltxt), path(strelka_snvs_fulltxt)

     output:
     tuple val(meta),
     path("${meta.lib}.unionSomaticVarsFull.txt")

     stub:
     """
       touch "${meta.lib}.unionSomaticVarsFull.txt"
     """
     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
     UnionSomaticCalls.pl ${strelka_snvs_fulltxt} ${strelka_indels_fulltxt} ${mutect_fulltxt} >${prefix}.unionSomaticVarsFull.txt
    
     """


}