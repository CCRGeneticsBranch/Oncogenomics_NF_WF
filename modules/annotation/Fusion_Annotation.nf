process Fusion_Annotation {

    tag "$meta.lib"
    //publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.id}/db", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),path(isoform_file),path(actionable_fusion),
    path(pfamdb),
    path(genome),
    val(genome_version_fusion_annotation),
    val(genome_version)

    output:
    tuple val(meta),
        path("${meta.lib}.annotated.txt")

    script:
     """
        fusionTools.py -i ${actionable_fusion} -m ${isoform_file} -o ${meta.lib}.annotated -t ${task.cpus} -p ${pfamdb} -f ${genome} -g /apps/data/gencode.${genome_version_fusion_annotation}.annotation.sorted.gtf.gz -n /apps/data/gencode.${genome_version_fusion_annotation}.canonical.txt -d /apps/data/gencode.${genome_version_fusion_annotation}.domains.tsv
     """
}


process Merge_fusion_annotation {

 tag "$meta.id"
 publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.id}/db", mode: "${params.publishDirMode}"

input:
 tuple val(meta),path(files),val(genome_version)

output:
 tuple val(meta),
 path("${meta.id}.fusion.txt"),
 path("${meta.id}.html")

stub:
 """
 touch "${meta.id}.fusion.txt",
 touch "${meta.id}.html"
 """

 script:

 """
 cat ${files[0]} |head -n 1 > ${meta.id}.fusion.txt

 for file in \${libs[@]}; do
    tail -n +2 \$file >> ${meta.id}.fusion.txt
 done

 python /apps/makeOutputHTML.py -i ${meta.id}.fusion.txt -o ${meta.id}.html -t /apps/data/template.html -c /apps/data/hg${genome_version}_cytoBand.txt
 """

}
