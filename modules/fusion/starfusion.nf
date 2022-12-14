process Starfusion{
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/starfusion", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        path(chimeric_junctions),
        path(genome_lib_dir)
    
    output:
    tuple val("${dataset_id}"),
        path("${dataset_id}.fusion_predictions.tsv")

    stub:
    """
    touch "${dataset_id}.fusion_predictions.tsv"
    """

    shell:
    '''
    set -exo pipefail
	STAR-Fusion \
        --genome_lib_dir !{genome_lib_dir} \
        --chimeric_junction !{chimeric_junctions} \
        --CPU !{task.cpus}
    cp STAR-Fusion_outdir/star-fusion.fusion_predictions.tsv !{dataset_id}.fusion_predictions.tsv
    '''
}