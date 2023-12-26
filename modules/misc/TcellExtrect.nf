process TcellExtrect {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/TCellExTRECT", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
     path(bam),
     path(index),
     path(bed),
     val(genome_version_tcellextrect)

     output:
     tuple val(meta),path("${meta.lib}_TCellExTRECT_naive.txt")

    stub:
    """
        touch "${meta.lib}_TCellExTRECT_naive.txt"
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"

    """

    Rscript /apps/runTCellExTRECT.R ${bam} ${bed} ${prefix} ${prefix}_TCellExTRECT ${genome_version_tcellextrect}
    """

}

process TcellExtrect_TN {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/TCellExTRECT", mode: "${params.publishDirMode}", pattern: "${meta.lib}*"
     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.normal_id}/TCellExTRECT", mode: "${params.publishDirMode}", pattern: "${meta.normal_id}*"
 

     input:
     tuple val(meta),
     path(Nbam),
     path(Nindex),
     path(Tbam),
     path(Tindex),
     path(bed),
     path(sequenza_alternate_sols),
     val(genome_version_tcellextrect)

     output:
     tuple val(meta),path("${meta.lib}*.txt"),path("${meta.normal_id}*.txt") , emit: naive_txt
     tuple val(meta),path("${meta.lib}*pdf"),path("${meta.normal_id}*.pdf") , emit: pdf , optional: true

    stub:
    """
        touch "${meta.lib}_TCellExTRECT_naive.txt"
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    def sequenza = sequenza_alternate_sols ? "${sequenza_alternate_sols}" : ""
    """
    purity=`head -2 ${sequenza_alternate_sols} | tail -1 | cut -f1`
    ploidy=`head -2 ${sequenza_alternate_sols} | tail -1 | cut -f2`
    Rscript /apps/runTCellExTRECT.R ${Tbam} ${bed} ${prefix} ${prefix}_TCellExTRECT ${genome_version_tcellextrect} \$purity \$ploidy
    Rscript /apps/runTCellExTRECT.R ${Nbam} ${bed} ${meta.normal_id} ${meta.normal_id}_TCellExTRECT ${genome_version_tcellextrect} \$purity \$ploidy
    """

}
