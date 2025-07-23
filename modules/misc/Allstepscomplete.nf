process Allstepscomplete {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}", mode: "${params.publishDirMode}"

    cache 'lenient'

    input:
    tuple val(meta),path(config)
    path ch_allcomplete

    output:

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
