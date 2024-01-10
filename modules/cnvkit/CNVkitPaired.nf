process CNVkitPaired {

 tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/cnvkit", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

    input:
    tuple val(meta),
    path(Nbam),
    path(Nindex),
    path(Tbam),
    path(Tindex),
    path(Tbed),
    path(cnv_ref_access),
    path(genome),
    path(genome_fai),
    path(genome_dict)

    output:
    tuple val(meta),path("${meta.lib}.cns"), emit: cnvkit_cns
    tuple val(meta),path("${meta.lib}.cnr"), emit: cnvkit_cnr
    tuple val(meta),path("${meta.lib}.pdf"), emit: cnvkit_pdf
    tuple val(meta),path("${meta.lib}_genelevel.txt"), emit: cnvkit_genelevel
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
    cnvkit.py batch -p ${task.cpus} --access ${cnv_ref_access} --fasta ${genome} --targets ${Tbed} ${Tbam} --output-dir . --normal ${Nbam}
    mv ${prefix}.final.cns ${prefix}.cns
    mv ${prefix}.final.cnr ${prefix}.cnr
    cnvkit.py scatter -s ${prefix}.cn{s,r} -o ${prefix}.pdf
    grep -v NOTFOUND NCI0439_T1D_E.cnr |grep -v Antitarget|perl -nae '\$F[3]=~s/_{2}.*//;print join("\t",@F)."\n"' > ${prefix}.cnr_filtered
    cnvkit.py genemetrics ${prefix}.cnr_filtered -t 0 -o ${prefix}_genelevel.txt

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    cnvkit: \$(cnvkit.py batch 2>&1  |grep -E '^CNVkit'|sed 's/CNVkit //')
END_VERSIONS


    """
}
