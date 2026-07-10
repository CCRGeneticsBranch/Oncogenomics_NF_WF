process Allstepscomplete {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}", mode: "${params.publishDirMode}"

    cache 'lenient'

    input:
    tuple val(meta),path(config)
    path complete_list

    output:

    path("successful*")

    stub:
    """
    touch "successful*"
    """

    script:
     def successful_file = params.genome_v == 'hg19' ? 'successful.txt' : 'successful_hg38.txt'
    """
     touch "${successful_file}"
    """

}
