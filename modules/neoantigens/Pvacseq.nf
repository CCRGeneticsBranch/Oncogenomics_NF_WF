process Split_vcf {

   tag "$meta.lib"

   publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/NeoAntigen", mode: "${params.publishDirMode}"

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


