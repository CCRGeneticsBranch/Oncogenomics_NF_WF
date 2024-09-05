process Mantis_MSI {

 tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

    input:
    tuple val(meta),
    path(Nbam),
    path(Nindex),
    path(Tbam),
    path(Tindex),
    path(loci_bed),
    path(genome)

    output:
    tuple val(meta),path("${meta.lib}_msi.txt"), emit: msi_txt
    tuple val(meta),path("${meta.lib}_msi.txt.status"), emit: msi_status

    stub:
    """
    touch "${meta.lib}_msi.txt"
    touch "${meta.lib}_msi.txt.status"

    """

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    mantis.py --bedfile ${loci_bed} \
            --genome  ${genome} \
            -n ${Nbam} -t ${Tbam} \
            -mrq 20.0 -mlq 25.0 \
            -mlc 20 -mrr 1 \
            -o ${meta.lib}_msi.txt \
            --threads ${task.cpus}

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    Mantis: \$(mantis.py 2>&1  |grep -E '^Microsatellite'|sed -e 's/.*(\\(.*\\)).*/\\1/')
END_VERSIONS


    """
}
