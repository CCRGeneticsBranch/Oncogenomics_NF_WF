process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_low'
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/qc", mode: "${params.publishDirMode}",pattern: "*config*.txt"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}",mode: "${params.publishDirMode}",pattern: "successful.txt"
    input:
    tuple val(meta), path(versions),val(pipeline_version)

    output:
    tuple val(meta), path("*config*.txt")    , emit: config
    //path "successful.txt"   , emit: successful

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """

    dumpsoftwareversions.py ${versions} ${nextflow.version} Software_versions ${meta.id} ${pipeline_version}
    sed -i "s/'//g" ${meta.id}.config*txt
    #touch successful.txt
    """
}
