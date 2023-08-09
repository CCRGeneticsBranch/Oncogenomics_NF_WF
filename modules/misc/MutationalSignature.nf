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