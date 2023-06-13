process Allstepscomplete {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}", mode: "${params.publishDirMode}"

    cache 'lenient'

    input:
    tuple val(meta),
        path(fusion),
        path(annotation),
        path(multiqc)

    output:
    tuple val(meta),
        path("successful.txt")

    stub:
    """
    touch "successful.txt"
    """

    script:
    """
     touch "successful.txt"

    """
}


