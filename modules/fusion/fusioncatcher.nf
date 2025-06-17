


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
    touch \$TMP/reads_ids_unmapped.txt
        fusioncatcher.py \
            -p ${task.cpus} \
            -d ${db} \
            -i ${trim[0]},${trim[1]} \
            -o \$TMP  2>&1 | tee fc_log.txt

        if grep -q "KeyError: 'Gene_1_symbol(5end_fusion_partner)'" fc_log.txt || grep -q "Auto-detect found SOLEXA FASTQ format" fc_log.txt; then
            echo "⚠️ No fusions detected. Creating dummy output file."
            touch ${prefix}.fusion-catcher.txt
            touch ${prefix}.summary_candidate_fusions.txt
        else
            cp \$TMP/final-list_candidate-fusion-genes.hg19.txt ${prefix}.fusion-catcher.txt
            cp \$TMP/summary_candidate_fusions.txt ${prefix}.summary_candidate_fusions.txt
        fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Fusioncatcher: \$(fusioncatcher --version|sed 's/.*fusioncatcher //')
    END_VERSIONS

    """

}
