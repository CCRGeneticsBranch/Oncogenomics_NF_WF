process OPENCRAVAT {
    tag "$meta.id"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/annotation/opencravat", mode: "${params.publishDirMode}"

    input:
    tuple val(meta), path(vcfs)

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

    /data/khanlab/projects/ngs_pipeline_testing/opencravat/oc_vg_2.13_install/bin/oc run \\
        ${vcfs.join(' ')} \\
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

process ADD_OPENCRAVAT_ANNOTATIONS_TN {
    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

    input:
    tuple val(meta), path(snpeff), path(filtered_tsv)

    output:
    tuple val(meta), path("${meta.lib}.HC_${meta.type}.annotated.txt")                        , emit: Tumor_hc_anno_txt
    tuple val(meta), path("${meta.normal_id}.HC_${meta.normal_type}.annotated.txt")           , emit: Normal_hc_anno_txt  , optional: true
    tuple val(meta), path("${meta.rna_lib}.HC_${meta.rna_type}.annotated.txt")                , emit: RNA_hc_anno_txt     , optional: true
    path "versions.yml"                                                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    if [[ ( "${meta.normal_type}" == "normal_DNA" || "${meta.normal_type}" == "blood_DNA" ) && "${meta.type}" == "tumor_DNA" && "${meta.rna_type}" != "tumor_RNA" ]]; then
        addOpenCravatAnnotations.py ${filtered_tsv} ${snpeff[0]} > ${meta.normal_id}.HC_${meta.normal_type}.annotated.txt
        addOpenCravatAnnotations.py ${filtered_tsv} ${snpeff[1]} > ${meta.lib}.HC_${meta.type}.annotated.txt
    elif [[ ( "${meta.normal_type}" == "normal_DNA" || "${meta.normal_type}" == "blood_DNA" ) && "${meta.type}" == "tumor_DNA" && "${meta.rna_type}" == "tumor_RNA" ]]; then
        addOpenCravatAnnotations.py ${filtered_tsv} ${snpeff[0]} > ${meta.normal_id}.HC_${meta.normal_type}.annotated.txt
        addOpenCravatAnnotations.py ${filtered_tsv} ${snpeff[1]} > ${meta.lib}.HC_${meta.type}.annotated.txt
        addOpenCravatAnnotations.py ${filtered_tsv} ${snpeff[2]} > ${meta.rna_lib}.HC_${meta.rna_type}.annotated.txt
    elif [[ ( "${meta.type}" == "tumor_DNA" || "${meta.type}" == "cell_line_DNA" ) && ( "${meta.rna_type}" == "tumor_RNA" || "${meta.rna_type}" == "cell_line_RNA" ) && ( "${meta.normal_type}" != "normal_DNA" && "${meta.normal_type}" != "blood_DNA" ) ]]; then
        addOpenCravatAnnotations.py ${filtered_tsv} ${snpeff[0]} > ${meta.lib}.HC_${meta.type}.annotated.txt
        addOpenCravatAnnotations.py ${filtered_tsv} ${snpeff[1]} > ${meta.rna_lib}.HC_${meta.rna_type}.annotated.txt
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
    END_VERSIONS
    """
}

process ADD_OPENCRAVAT_ANNOTATIONS_SOMATIC_VARIANTS {
    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
    path(mutect_txt),
    path(strelka_indels_txt),
    path(strelka_snvs_txt),
    path(filtered_tsv)

    output:
    tuple val(meta),
    path("${meta.lib}.MuTect.annotated.txt"),
    path("${meta.lib}.strelka.indels.annotated.txt"),
    path("${meta.lib}.strelka.snvs.annotated.txt")

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    addOpenCravatAnnotations.py ${filtered_tsv} ${mutect_txt}         > ${prefix}.MuTect.annotated.txt
    addOpenCravatAnnotations.py ${filtered_tsv} ${strelka_indels_txt} > ${prefix}.strelka.indels.annotated.txt
    addOpenCravatAnnotations.py ${filtered_tsv} ${strelka_snvs_txt}   > ${prefix}.strelka.snvs.annotated.txt

    """
}
