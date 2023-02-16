process MergeHLA {
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/HLA", mode: "${params.publishDirMode}"
    input:
//    path(mergeHLA_script)
    tuple val(dataset_id),
        path(seq2hla),
        path(hlaminer)

    output:
	path "${dataset_id}.Calls.txt"

    stub:
    """
    touch "${dataset_id}.Calls.txt"
    """


    shell:
    '''
    set -exo pipefail

    consensusHLA.pl !{hlaminer} !{seq2hla} |sort > !{dataset_id}.Calls.txt

    '''
}


