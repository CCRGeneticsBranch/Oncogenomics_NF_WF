process Optitype {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/HLA", mode: "${params.publishDirMode}", pattern: "${meta.lib}*"

    input:
    tuple val(meta), path(trim)

    output:
    tuple val(meta), path("${meta.lib}_optitype.txt"), emit: Optitype_output
    path "versions.yml"             , emit: versions

    stub:
    """
    touch "${meta.lib}_optitype.txt"
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    if [[ "${meta.type}" == "tumor_DNA" || "${meta.type}" == "cell_line_DNA" || "${meta.type}" == "normal_DNA" || "${meta.type}" == "blood_DNA" ]]; then
        type="--dna"
    elif [[ "${meta.type}" == "tumor_RNA" || "${meta.type}" == "cell_line_RNA" ]]; then
        type="--rna"
    fi
    TMP=tmp
    mkdir \$TMP
    trap 'rm -rf "\$TMP"' EXIT
    if [[ \$(OptiTypePipeline.py -i ${trim[0]} ${trim[1]} \$type -v -o \$TMP 2>&1 | grep -m 1 -e 'need more than 0 values to unpack' -e 'The constraint expression resolved to a trivial Boolean (False) instead of a Pyomo object.') ]]; then
        echo "no HLA called, output empty file"
        touch "${meta.lib}_optitype.txt"
    else
        cp \$TMP/*/*_result.tsv ${meta.lib}_optitype.txt
    fi
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Optitype: \$(echo \$SINGULARITY_NAME | sed 's/fred2-optitype-release-//g'|sed 's/.img//g')
    END_VERSIONS
    """
}


process HLA_HD {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/HLA", mode: "${params.publishDirMode}", pattern: "${meta.lib}*"

    input:
    tuple val(meta), path(trim)

    output:
    tuple val(meta), path("${meta.lib}_HLAHD.txt"), emit: hlahd_output
    path "versions.yml"             , emit: versions

    stub:
    """
    touch "${meta.lib}_HLAHD.txt"
    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    TMP=tmp
    mkdir \$TMP
    trap 'rm -rf "\$TMP"' EXIT
    gunzip -c ${trim[0]} > ${meta.lib}_R1_unzipped.fastq
    gunzip -c ${trim[1]} > ${meta.lib}_R2_unzipped.fastq

    hlahd.sh -t ${task.cpus} \
        -m 70 -f /data2/hlahd.1.7.0/freq_data \
        ${meta.lib}_R1_unzipped.fastq ${meta.lib}_R2_unzipped.fastq \
        /data2/hlahd.1.7.0/HLA_gene.split.txt \
        /data2/hlahd.1.7.0/dictionary \
        ${meta.lib} \
        \$TMP
    cp \$TMP/${meta.lib}/result/${meta.lib}_final.result.txt ${meta.lib}_HLAHD.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        HLA-HD: \$(hlahd.sh 2>&1 | grep -o 'HLA-HD version [0-9.]*' | awk '{print \$NF}')
    END_VERSIONS
    """
}

process Merge_new_HLA {

    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/HLA", mode: "${params.publishDirMode}"
    input:

    tuple val(meta),
        path(optitype),
        path(hlahd)

    output:
	tuple val(meta),path("${meta.lib}.Calls.txt")

    stub:
    """
    touch "${meta.lib}.Calls.txt"
    """


    script:
    """
    HLA_consensus.py ${optitype} ${hlahd} ${meta.lib}.Calls.txt

    """
}
