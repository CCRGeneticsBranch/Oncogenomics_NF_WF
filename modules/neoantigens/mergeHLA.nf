process MergeHLA {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/HLA", mode: "${params.publishDirMode}"
    input:

    tuple val(meta),
        path(seq2hla),
        path(hlaminer)

    output:
	tuple val(meta),path("${meta.lib}.Calls.txt")

    stub:
    """
    touch "${meta.lib}.Calls.txt"
    """


    script:
    """
    set -exo pipefail

    consensusHLA.pl  ${seq2hla} ${hlaminer} |sort > ${meta.lib}.Calls.txt

    """
}


