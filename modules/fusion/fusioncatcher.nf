


process Fusioncatcher {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/fusion", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

    input:
    tuple val(meta),path(trim),path(db)

    output:
    tuple val(meta),path("${meta.lib}.fusion-catcher.txt"),  emit : fc_output
    tuple val(meta),path("${meta.lib}.summary_candidate_fusions.txt"), emit : fc_summary
    path "versions.yml"             , emit: versions

    stub:
    """
    touch "${meta.lib}.fusion-catcher.txt"
    touch "${meta.lib}.summary_candidate_fusions.txt"
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    TMP=tmp/
    mkdir -p \$TMP
    trap 'rm -rf "\$TMP"' EXIT

        fusioncatcher.py \
            -p ${task.cpus} \
            -d ${db} \
            -i ${trim[0]},${trim[1]} \
            -o \$TMP

        cp \$TMP/final-list_candidate-fusion-genes.hg19.txt ${prefix}.fusion-catcher.txt
        cp \$TMP/summary_candidate_fusions.txt ${prefix}.summary_candidate_fusions.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Fusioncatcher: \$(fusioncatcher --version|sed 's/.*fusioncatcher //')
    END_VERSIONS

    """

}
