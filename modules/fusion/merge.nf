process Mergefusion {
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/Actionable", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        val(library),
        path(arriba), 
        path(FC),
        path(FC_summary),
        path(SF)
    
    output:
    tuple val("${dataset_id}"),
        val("${library}"),
        path("${dataset_id}.fusion.actionable.txt")

    stub:
    """
    touch "${dataset_id}.fusion.actionable.txt"
    """

    shell:
    '''
    ActionableFusion.v1.pl  !{library} !{FC} !{SF} !{arriba} $PWD | awk 'NR<2{print $0;next}{print $0| "sort "}' > !{dataset_id}.fusion.actionable.txt

    '''
}
