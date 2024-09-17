process Combine_variants  {

   tag "$meta.lib"

   publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/NeoAntigen", mode: "${params.publishDirMode}", pattern: "*.txt"

  input:

  tuple val(meta),path(mutect_raw_vcf),path(strelka_indel_raw_vcf),path(strelka_snvs_raw_vcf)

  output:
  tuple val(meta),path("${meta.lib}.final.vcf.tmp") , emit: combined_vcf_tmp

  stub:
  """
  touch "${meta.lib}.final.vcf.tmp"

  """

  script:
  def prefix = task.ext.prefix ?: "${meta.lib}"
  """
  consensusSomaticVCF.pl -vcf ${strelka_indel_raw_vcf},${strelka_snvs_raw_vcf},${mutect_raw_vcf} -order ${meta.normal_id},${prefix} -filter REJECT |vcf-subset -u -c ${prefix} > ${prefix}.final.vcf.tmp

  """

}

process VEP {

   tag "$meta.lib"

   publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/NeoAntigen", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

  input:
  tuple val(meta),path(combined_vcf_tmp),path(vep_cache)

  output:
  tuple val(meta),path("${meta.lib}_final.vcf") , emit: vep_out
  path "versions.yml"             , emit: versions

  stub:
  """
  touch "${meta.lib}_final.vcf"
  """

  script:
  def prefix = task.ext.prefix ?: "${meta.lib}"
  """
  /opt/vep/src/ensembl-vep/vep -i ${combined_vcf_tmp} --format vcf --plugin Downstream --plugin Wildtype \
                  --terms SO --offline --cache --dir ${vep_cache} \
                   --assembly GRCh37 \
                  --output_file ${prefix}_final.vcf --vcf --force_overwrite --no_check_variants_order

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      VEP: \$(/opt/vep/src/ensembl-vep/vep|grep ensembl-vep|sed 's/ensembl-vep          : //')
  END_VERSIONS
  """



}
