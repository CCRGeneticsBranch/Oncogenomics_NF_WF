process Optitype {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/HLA", mode: "${params.publishDirMode}"

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
    if [[ "${meta.type}" == "tumor_DNA" || "${meta.type}" == "cell_line_DNA" || "${meta.type}" == "normal_DNA" ]]; then
        type="--dna"
    elif [[ "${meta.type}" == "tumor_RNA" || "${meta.type}" == "cell_line_RNA" ]]; then
        type="--rna"
    fi
    [ -d "output_optitype" ] && rm -rf "output_optitype" ; mkdir "output_optitype"
    OptiTypePipeline.py -i ${trim[0]} ${trim[1]} \$type -v -o output_optitype
    cp output_optitype/*/*_result.tsv ${meta.lib}_optitype.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Optitype: \$(echo \$SINGULARITY_NAME | sed 's/optitype_release-//g'|sed 's/.sif//g')
    END_VERSIONS
    """
}


process HLA_HD {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/HLA", mode: "${params.publishDirMode}"

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
    gunzip -c ${trim[0]} > ${meta.lib}_R1_unzipped.fastq
    gunzip -c ${trim[1]} > ${meta.lib}_R2_unzipped.fastq
    [ -d "output_hlahd" ] && rm -rf "output_hlahd" ; mkdir "output_hlahd"
    hlahd.sh -t ${task.cpus} \
        -m 70 -f /data2/hlahd.1.7.0/freq_data \
        ${meta.lib}_R1_unzipped.fastq ${meta.lib}_R2_unzipped.fastq \
        /data2/hlahd.1.7.0/HLA_gene.split.txt \
        /data2/hlahd.1.7.0/dictionary \
        ${meta.lib} \
        output_hlahd
    cp output_hlahd/${meta.lib}/result/${meta.lib}_final.result.txt ${meta.lib}_HLAHD.txt

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
	tuple val(meta),path("${meta.lib}.newCalls.txt")

    stub:
    """
    touch "${meta.lib}.newCalls.txt"
    """


    script:
    """
    HLA_consensus.py ${optitype} ${hlahd} ${meta.lib}.newCalls.txt

    """
}
