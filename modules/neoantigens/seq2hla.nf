process Seq2HLA {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/HLA/seq2HLA", mode: "${params.publishDirMode}"

    input:
    tuple val(meta), path(trim)

    output:
    tuple val(meta), path("${meta.lib}-ClassI.HLAgenotype4digits")

    stub:
    """
    touch "${meta.lib}-ClassI.HLAgenotype4digits"
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    set -exo pipefail
    seq2HLA \
        references/ \
        -1 ${trim[0]} \
        -2 ${trim[1]} \
        -p ${task.cpus} \
        -r "${prefix}"
    """
}

process Seq2HLA_exome {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/HLA/seq2HLA", mode: "${params.publishDirMode}"

    input:
    tuple val(meta), path(r1fq),path(r2fq)

    output:
    tuple val(meta), path("${meta.lib}-ClassI.HLAgenotype4digits"), emit: seq2hla_output
    path "versions.yml"             , emit: versions

    stub:
    """
    touch "${meta.lib}-ClassI.HLAgenotype4digits"
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    set -exo pipefail
    seq2HLA \
        references/ \
        -1 ${r1fq} \
        -2 ${r2fq} \
        -p ${task.cpus} \
        -r "${prefix}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Seq2HLA: \$(seq2HLA --version |sed 's/.*seq2HLA.py //')
    END_VERSIONS
    """
}
