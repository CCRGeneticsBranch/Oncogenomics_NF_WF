process Seq2HLA {
    tag { dataset_id }

    publishDir "$params.resultsdir/$dataset_id/Seq2HLA", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        path(r1), 
        path(r2)
    
    output:
    tuple val("${dataset_id}"),
        path("${dataset_id}-ClassI.HLAgenotype4digits")

    stub:
    """
    touch "${dataset_id}-ClassI.HLAgenotype4digits"
    """

    shell:
    '''
    set -exo pipefail

    seq2HLA \
        references/ \
        -1 !{r1} \
        -2 !{r2} \
        -p !{task.cpus} \
        -r "!{dataset_id}"
    '''
}


