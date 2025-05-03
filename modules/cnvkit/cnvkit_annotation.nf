process CNVkitAnnotation {

    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/cnvkit", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"
    input:
    tuple val(meta),path(target_bed),path(cnvkit_call_cns)
    path(combined_gene_list)
    output:
    tuple val(meta),path("${meta.lib}_genelevel.txt"), emit: cnvkit_genelevel

    stub:
    """
    touch "${meta.lib}.genelevel.txt"
    """
    script:
    """
    awk -F'\t' 'NF==3 {\$4="NOTFOUND___0"; \$5="NULL"} {OFS="\t"; print}' ${target_bed} > 5col_target.bed
    cnvkit_geneann_ngs.py  ${meta.lib} ./ 5col_target.bed  ${combined_gene_list} ${cnvkit_call_cns}
    """
}
