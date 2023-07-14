process Fusion_Annotation {

    tag "$meta.id"
//    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.id}/db", mode: "${params.publishDirMode}"

    input:
    path(isoform_files)
    tuple val(meta), path(fusion_actionable),path(pfamdb),path(genome)
    val(genome_version_fusion_annotation)
    val(genome_version)

    output:
    tuple val(meta),
        path("${meta.id}.check")


    script:
     """
x=${isoform_files.join(',')}
IFS=',' read -ra items <<< \$x
for item in \${items[@]}; do
  echo "Item: \$item"
  fusionTools.py -i ${fusion_actionable} -m \$item -o \${item}_fusion.txt -t ${task.cpus} -p ${pfamdb} -f ${genome} -g /apps/data/gencode.${genome_version_fusion_annotation}.annotation.sorted.gtf.gz -n /apps/data/gencode.${genome_version_fusion_annotation}.canonical.txt -d /apps/data/gencode.${genome_version_fusion_annotation}.domains.tsv
done
cat *fusion.txt > ${meta.id}.fusion.txt
/apps/makeOutputHTML.py -i ${meta.id}.fusion.txt -o ${out_file}.html -t /opt/data/template.html -c /apps/data/hg${genome_version}_cytoBand.txt

     """

}