process Combine_variants  {

   tag "$meta.lib"

   publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/NeoAntigen", mode: "${params.publishDirMode}", pattern: "*.txt"

  input:

  tuple val(meta),path(mutect_raw_vcf),path(strelka_indel_raw_vcf),path(strelka_snvs_raw_vcf)
  tuple val(meta2),path(hlaminer),path(seq2hla)
  
  output:
  tuple val(meta),path("${meta.lib}.final.vcf.tmp")
  tuple val(meta2),path("${meta2.lib}.Calls.txt")

  stub:
  """
  touch "${meta.lib}.final.vcf.tmp"
  touch "${meta2.lib}.Calls.txt"
  """

  script:
  def prefix = task.ext.prefix ?: "${meta.lib}"
  """
  consensusSomaticVCF.pl -vcf ${strelka_indel_raw_vcf},${strelka_snvs_raw_vcf},${mutect_raw_vcf} -order ${meta2.lib},${prefix} -filter REJECT |vcf-subset -u -c ${prefix} > ${prefix}.final.vcf.tmp
  consensusHLA.pl ${seq2hla} ${hlaminer} |sort > ${meta2.lib}.Calls.txt

  """

}