process MutationBurden {

    tag "$meta.id"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),path(Mutect_annotationfull),path(strelka_indels_annotationfull),path(strelka_snvs_annotationfull)
    path(tumor_target_capture)
    tuple val(meta),val(normal),val(tumor),val(vaf)
    val(mutect_ch)
    val(strelka_indelch)
    val(strelka_snvsch)

    output:
    tuple val(meta),
    path("${meta.lib}.*mutationburden.txt")

    script:
    """
    mutationBurden.pl  ${strelka_indels_annotationfull} ${tumor_target_capture} ${tumor} ${normal} ${vaf} > ${meta.lib}.${strelka_indelch}.mutationburden.txt
    
    """

}