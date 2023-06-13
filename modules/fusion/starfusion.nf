process Starfusion{
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/fusion", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),path(chimeric_junctions),path(genome_lib_dir)

    output:
    tuple val(meta),path("${meta.lib}.STAR-fusion.txt")

    stub:
    """
    touch ("${meta.lib}.STAR-fusion.txt")
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    set -exo pipefail
        STAR-Fusion \
        --genome_lib_dir ${genome_lib_dir} \
        --chimeric_junction ${chimeric_junctions} \
        --CPU ${task.cpus}
    cp STAR-Fusion_outdir/star-fusion.fusion_predictions.tsv ${prefix}.STAR-fusion.txt
    """
}
