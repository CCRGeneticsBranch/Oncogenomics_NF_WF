process HLAminer {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/HLA/HLAminer", mode: "${params.publishDirMode}"

    input:
    tuple val(meta), path(trim)

    output:
    tuple val(meta), path("${meta.lib}_HLAminer_HPTASR.csv")


    stub:
    """
    touch "${meta.lib}_HLAminer_HPTASR.csv"
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    set -exo pipefail

    echo ${trim[0]} >patient.fof

    echo ${trim[1]} >>patient.fof

    HPTASRrnaseq_classI.sh .
    mv HLAminer_HPTASR.csv ${prefix}_HLAminer_HPTASR.csv
    """
}

process HLAminer_exome {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/HLA/HLAminer", mode: "${params.publishDirMode}"

    input:
    tuple val(meta), path(r1fq),path(r2fq)

    output:
    tuple val(meta), path("${meta.lib}_HLAminer_HPTASR.csv"), emit: hlaminer_output
    path "versions.yml"             , emit: versions

    stub:
    """
    touch "${meta.lib}_HLAminer_HPTASR.csv"
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    set -exo pipefail

    echo ${r1fq} >patient.fof

    echo ${r2fq} >>patient.fof

    HPTASRwgs_classI.sh .
    mv HLAminer_HPTASR.csv ${prefix}_HLAminer_HPTASR.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        HLAminer: \$(HPTASRwgs_classI.sh |awk -F '[][]' '/^Running:/{print \$2}')
    END_VERSIONS
    """
}
