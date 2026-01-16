process OPENCRAVAT {
    tag "$meta.id"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/annotation/opencravat", mode: "${params.publishDirMode}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.tsv")     , emit: tsv
    tuple val(meta), path("*.html")    , emit: html, optional: true
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """

    #/data/khanlab/projects/ngs_pipeline_testing/opencravat/oc_vg_2.13_install/bin/oc run \\
    oc run \\
        $vcf \\
        -a gnomad3 \\
        -l ${params.genome_v} \\
        -n ${prefix} \\
        -t text \\
        --mp ${task.cpus} \\
        --md ${params.modules_oc} \\
        -d . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        opencravat: \$(oc version 2>&1 | grep -oP 'Open-CRAVAT \\K[0-9.]+')
    END_VERSIONS
    """
}

process FILTER_OPENCRAVAT_OUTPUT {
    tag "$meta.id"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/annotation/opencravat", mode: "${params.publishDirMode}"

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("${meta.id}.Annotations.coding.rare.txt"), emit: rare
    tuple val(meta), path("${meta.id}.Annotations.final.txt"), emit: full

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    filter_opencravat_output.py \\
        ${tsv} \\
        --sample-id ${meta.id} \\
        $args

    """
}

process ADD_OPENCRAVAT_ANNOTATIONS {
    tag "$meta.id"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/annotation/opencravat", mode: "${params.publishDirMode}"

    input:
    tuple val(meta), path(snpeff_txt),path(filtered_tsv)

    output:
    tuple val(meta), path("*.annotated.txt"), emit: annotated_txt


    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    addOpenCravatAnnotations.py \\
        ${filtered_tsv} \\
        ${snpeff_txt} \\
        > ${meta.lib}.HC_${meta.type}.annotated.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}
