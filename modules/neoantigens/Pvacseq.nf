process Split_vcf {

   tag "$meta.lib"

   //publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/NeoAntigen", mode: "${params.publishDirMode}"

  input:
  tuple val(meta),path(vep_vcf)

  output:
  tuple val(meta),path("split/*vcf")

  stub:
  """
  touch "split/*vcf"
  """

  script:
  def prefix = task.ext.prefix ?: "${meta.lib}"
  """
  mkdir -p split
  split_vcf.py ${vep_vcf} split/
  """

}


process Pvacseq {
   tag "$meta.lib"

   //publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/NeoAntigen/split", mode: "${params.publishDirMode}"

  input:
  tuple val(meta),path(split_vcf),path(combined_HLA_normalsamples)

  output:
  tuple val(meta),path("*filtered.tsv"), emit: pvacseq_output_ch
  path "versions.yml"             , emit: versions

  stub:
  """
  touch "*filtered.tsv"
  """
  script:
  def prefix = task.ext.prefix ?: "${meta.lib}"
  """
  SAMPLE=`basename ${split_vcf} .vcf`
  mkdir \$SAMPLE
  pvacseq_script.sh ${combined_HLA_normalsamples} ${split_vcf} \$SAMPLE \$SAMPLE ${task.cpus}

  if [ -f \$SAMPLE/MHC_Class_I/*filtered.tsv ]
  then
    cp \$SAMPLE/MHC_Class_I/*filtered.tsv .
  else
    touch \$SAMPLE.filtered.tsv
  fi

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    pvacseq: \$(pvactools --version)
END_VERSIONS

  """

}

process Merge_Pvacseq_vcf {
  tag "$meta.lib"

  publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/NeoAntigen", mode: "${params.publishDirMode}"

  input:
  tuple val(meta), path(vcf_files)

  output:
  tuple val(meta),path("${meta.lib}.final.txt")

  stub:
  """
  touch "${meta.lib}.final.txt"
  """

  script:
  def prefix = task.ext.prefix ?: "${meta.lib}"
  """
  awk 'FNR>1 || NR==1' ${vcf_files.join(' ')} > ${meta.lib}.filtered.tsv
  awk 'NR == 1; NR > 1 {print \$0 | "sort -n"}' ${meta.lib}.filtered.tsv|uniq > ${meta.lib}.final.uniq.tsv
  process_pVACSeq.pl ${meta.lib}.final.uniq.tsv |awk 'NR == 1; NR > 1 {print \$0 |"sort -n"}'|uniq > ${meta.lib}.final.txt
  """


}
