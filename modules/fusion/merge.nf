process Mergefusion {
    tag "$meta.lib"

    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/Actionable", mode: "${params.publishDirMode}"

    cache 'lenient'

    input:
    tuple val(meta),path(arriba),path(FC),path(SF)

    output:
    tuple val(meta),path("${meta.lib}.actionable.fusion.txt")

    stub:
    """
    touch "${meta.lib}.actionable.fusion.txt"
    """

    script:
    """
    ActionableFusion.v1.pl  ${meta.lib} ${FC} ${SF} ${arriba}  | awk 'NR<2{print \$0;next}{print \$0| "sort "}' > ${meta.lib}.actionable.fusion.txt

    """
}
