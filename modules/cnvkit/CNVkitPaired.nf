process CNVkitPaired {

 tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/cnvkit", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

    input:
    tuple val(meta),path(Nbam),path(Nindex),path(Tbam),path(Tindex),path(Tbed),path(sequenza_alternate),path(mutect_raw_vcf)
    path(cnv_ref_access)
    path(genome)
    path(genome_fai)
    path(genome_dict)

    output:
    tuple val(meta),path("${meta.lib}.cns"), emit: cnvkit_cns
    tuple val(meta),path("${meta.lib}.cnr"), emit: cnvkit_cnr
    tuple val(meta),path("${meta.lib}.pdf"), emit: cnvkit_pdf
    tuple val(meta),path("${meta.lib}.call.cns"), emit: cnvkit_call_cns
    tuple val(meta),path("${meta.lib}.call.cns_lowPurity"), emit: cnvkit_call_cns_lowPurity , optional: true
    //tuple val(meta),path("${meta.lib}_genelevel.txt"), emit: cnvkit_genelevel
    path "versions.yml"             , emit: versions

    stub:
     """
     touch "${meta.lib}.cns"
     touch "${meta.lib}.cnr"
     touch "${meta.lib}.pdf"
     touch "${meta.lib}_genelevel.txt"

     """
    script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    cnvkit.py batch -p ${task.cpus} --access ${cnv_ref_access} --fasta ${genome} --targets ${Tbed} ${Tbam} --drop-low-coverage --output-dir . --normal ${Nbam}
    mv ${prefix}.final.cns ${prefix}.cns
    mv ${prefix}.final.cnr ${prefix}.cnr
    cnvkit.py scatter -s ${prefix}.cn{s,r} -o ${prefix}.pdf
    if grep -q cellularity ${sequenza_alternate}; then
        PURITY=`awk '{ print \$1 }' ${sequenza_alternate} | sed -n '2p'`
        echo "purity is  \$PURITY"
        cnvkit.py call -m clonal ${prefix}.cns --purity \$PURITY -v ${mutect_raw_vcf} -o ${prefix}.call.cns
        limit=0.3
        #if [ 1 -eq "\$(echo "\$PURITY < \$limit" | bc)" ]; then
        if awk "BEGIN {exit !(\$PURITY < \$limit)}"; then
            mv ${prefix}.call.cns ${prefix}.call.cns_lowPurity
            touch ${prefix}.call.cns
        fi
    else
        touch ${prefix}.call.cns
    fi

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    cnvkit: \$(cnvkit.py batch 2>&1  |grep -E '^CNVkit'|sed 's/CNVkit //')
END_VERSIONS


    """
}
