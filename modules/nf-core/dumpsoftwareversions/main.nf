process CUSTOM_DUMPSOFTWAREVERSIONS {
    label 'process_low'

    input:
    path versions

    output:
    path "software_versions.yml"    , emit: yml
    path "software_versions_mqc.yml", emit: mqc_yml
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    cat ${versions}
    echo ${task.process}
    dumpsoftwareversions.py
    """
}
