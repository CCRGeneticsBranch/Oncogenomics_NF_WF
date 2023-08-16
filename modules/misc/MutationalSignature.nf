process MutationalSignature {

    tag "$meta.id"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/Actionable", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),path(unionSomaticVarsFull)

    output:
    tuple val(meta),
    path("${meta.lib}.mutationalSignature.pdf")

    stub:
    """   
    touch "${meta.lib}.mutationalSignature.pdf"
    """    
    
    script:
    """
    awk '{OFS="\t"}{print \$1,\$2,\$4,\$5,"${meta.lib}"}' ${unionSomaticVarsFull} |sed -e '1s/${meta.lib}/Sample/g' > ${meta.lib}.mutationalSignature.pdf.tmp
    mutationSignature.R --input ${meta.lib}.mutationalSignature.pdf.tmp --sample ${meta.lib} --output ${meta.lib}.mutationalSignature.pdf
    rm -rf ${meta.lib}.mutationalSignature.pdf.tmp
    """

}


process Cosmic3Signature {
    tag "$meta.id"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/Actionable/input", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),path(mutect),path(strelka_indels),path(strelka_snvs)
    path(cosmic_indel_rda)
    path(cosmic_genome_rda)
    path(cosmic_dbs_rda)

    output:
    tuple val(meta),
    path("${meta.lib}.Indel83_cosmic_v3.pdf"),
    path("${meta.lib}.SBS96_cosmic_v3.pdf"),
    path("${meta.lib}.DBS78_cosmic_v3.pdf")


    stub:
    """   
    touch "${meta.lib}.Indel83_cosmic_v3.pdf"
    touch "${meta.lib}.SBS96_cosmic_v3.pdf"
    touch "${meta.lib}.DBS78_cosmic_v3.pdf"
    """    

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    if [ -d "input"] ; then rm -rf input ; fi
    mkdir input
    mv ${mutect} ${strelka_indels} ${strelka_snvs} ./input
    matrixgenerator.py ${prefix} ./input
    cut -f1 ./input/output/ID/${prefix}.ID83.all |awk -F ":" -v OFS="_" '{print \$2,\$3,\$1,\$4}'|sed  s'/Del/DEL/g'|sed s'/Ins/INS/g' |sed s'/__MutationType_/MutationType/g'|sed s'/_R_/_repeats_/g'|sed s'/_M_/_MH_/g'| sed s'/5/5+/g' |paste - ./input/output/ID/${prefix}.ID83.all|cut -f1,3 > ./input/output/ID/${prefix}.ID83.all_updatedcolumns
    deconstructsigs_indels.R ./input/output/ID/${prefix}.ID83.all_updatedcolumns ${prefix}.Indel83 ${cosmic_indel_rda}
    deconstructsigs_sbs.R ./input/output/SBS/${prefix}.SBS96.all ${prefix}.SBS96 ${cosmic_genome_rda}
    deconstructsigs_dbs.R ./input/output/DBS/${prefix}.DBS78.all ${prefix}.DBS78 ${cosmic_dbs_rda}

    """ 

}

