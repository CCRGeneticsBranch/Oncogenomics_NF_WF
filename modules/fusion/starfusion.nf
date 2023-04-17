process Starfusion{
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/${library}/fusion", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        val(library),
        path(chimeric_junctions),
        path(genome_lib_dir)

    output:
    tuple val("${dataset_id}"),
        val("${library}"),
        path("${library}.STAR-fusion.txt")

    stub:
    """
    touch ("$library}.STAR-fusion.txt")
    """

    shell:
    '''
    set -exo pipefail
        STAR-Fusion \
        --genome_lib_dir !{genome_lib_dir} \
        --chimeric_junction !{chimeric_junctions} \
        --CPU !{task.cpus}
    cp STAR-Fusion_outdir/star-fusion.fusion_predictions.tsv !{library}.STAR-fusion.txt
    '''
}
