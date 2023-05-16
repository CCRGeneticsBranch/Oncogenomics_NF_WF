process Allstepscomplete {
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/${params.casename}", mode: "${params.publishDirMode}"

    cache 'lenient'

    input:
    tuple val(dataset_id),
        val(library),
        path(fusion),
        path(annotation)

    output:
    tuple val("${dataset_id}"),
        val("${library}"),
        path("successful.txt")

    stub:
    """
    touch "successful.txt"
    """

    shell:
    '''
     touch "successful.txt"

    '''
}


