process Mutect {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(Tbam),path(Tindex),path(Tbed)
     tuple val(meta2),path(Nbam),path(Nindex)
     path genome
     path genome_fai
     path genome_dict
     path dbsnp_138_b37_vcf
     path cosmic_v67_hg19_vcf


     output:
     tuple val(meta),path("${meta.lib}.MuTect_raw.vcf"),    emit: mutect_raw_vcf
     tuple val(meta), path("${meta.lib}.mutect.call_stats.txt"),          emit: mutect_stats
     tuple val(meta),path("${meta.lib}.mutect.coverage.wig.txt"),        emit: coverage_wig

     stub:
     """
        touch "${meta.lib}.MuTect_raw.vcf"
        touch "${meta.lib}.mutect.call_stats.txt"
        touch "${meta.lib}.mutect.coverage.wig.txt"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
     java -Xmx10g  -jar /opt/mutect-1.1.7.jar -T MuTect  --reference_sequence ${genome} \
            --cosmic ${cosmic_v67_hg19_vcf} \
            --dbsnp ${dbsnp_138_b37_vcf} \
            --input_file:normal ${Nbam} \
            --input_file:tumor ${Tbam} \
            --intervals ${Tbed} \
            --coverage_file ${prefix}.mutect.coverage.wig.txt \
            --out  ${prefix}.mutect.call_stats.txt \
            --vcf ${prefix}.MuTect_raw.vcf \
            --max_alt_allele_in_normal_fraction 0.05 \
            --max_alt_alleles_in_normal_count 4 --min_qscore 20 -rf MappingQuality -mmq 30 

     """
}


process Mutect_order {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),path(mutect_vcf)

    output:
     tuple val(meta),path("${meta.lib}.MuTect.raw.vcf")

    stub:
     """
     touch "${meta.lib}.MuTect.raw.vcf
     """
    
    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
     """

    vcfOrderCol.R -i ${mutect_vcf}  -o ${prefix}.MuTect.raw.vcf

     """

}

