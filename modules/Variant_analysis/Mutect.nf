process Mutect {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

     input:
     tuple val(meta),
     path(Nbam),
     path(Nindex),
     path(Tbam),
     path(Tindex),
     path(Tbed),
     path(genome),
     path(genome_fai),
     path(genome_dict),
     path(dbsnp_138_b37_vcf),
     path(cosmic_v67_hg19_vcf)


     output:
     tuple val(meta),path("${meta.lib}.MuTect_raw.vcf"),    emit: mutect_raw_vcf
     tuple val(meta), path("${meta.lib}.mutect.call_stats.txt"),          emit: mutect_stats
     tuple val(meta),path("${meta.lib}.mutect.coverage.wig.txt"),        emit: coverage_wig
     path "versions.yml"             , emit: versions

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

     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
         mutect: \$(java  -jar /opt/mutect-1.1.7.jar -T MuTect -version)
     END_VERSIONS
     """
}


process Mutect2 {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

     input:
     tuple val(meta),
     path(Nbam),
     path(Nindex),
     path(Tbam),
     path(Tindex),
     path(Tbed)
     path genome
     path genome_fai
     path genome_dict
     path germline_resource
     path germline_resource_idx
     path pon
     path pon_idx

     output:
     tuple val(meta),path("${meta.lib}.Mutect2_raw.vcf.gz"),          emit: mutect2_raw_vcf
     tuple val(meta),path("${meta.lib}.Mutect2_raw.vcf.gz.tbi"),      emit: mutect2_raw_vcf_tbi
     tuple val(meta),path("${meta.lib}.Mutect2_raw.vcf.gz.stats"),    emit: mutect2_stats
     tuple val(meta),path("${meta.lib}.f1r2.tar.gz"),                 emit: f1r2
     path "versions.yml"             , emit: versions

     stub:
     """
        touch "${meta.lib}.Mutect2_raw.vcf.gz"
        touch "${meta.lib}.Mutect2_raw.vcf.gz.tbi"
        touch "${meta.lib}.Mutect2_raw.vcf.gz.stats"
        touch "${meta.lib}.f1r2.tar.gz"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"

     """
     module load GATK/4.6.2.0
     gatk --java-options "-Xmx40g" Mutect2 \
            -R ${genome} \
            -I ${Tbam} \
            -I ${Nbam} \
            -normal ${meta.normal_id} \
            --germline-resource ${germline_resource} \
            --panel-of-normals ${pon} \
            -L ${Tbed} \
            --f1r2-tar-gz ${prefix}.f1r2.tar.gz \
            -O ${prefix}.Mutect2_raw.vcf.gz

     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
         gatk: \$(gatk --version 2>&1 | grep -oP 'v[0-9.]+')
     END_VERSIONS
     """
}

process learn_read_orientation_model {

     tag "$meta.lib"

     input:
     tuple val(meta), path(f1r2)

     output:
     tuple val(meta), path("${meta.lib}.read_orientation_model.tar.gz"), emit: orientation_model
     path "versions.yml",                                                 emit: versions

     stub:
     """
        touch "${meta.lib}.read_orientation_model.tar.gz"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
     module load GATK/4.6.2.0
     gatk --java-options "-Xmx16g" LearnReadOrientationModel \
          -I ${f1r2} \
          -O ${prefix}.read_orientation_model.tar.gz

     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
         gatk: \$(gatk --version 2>&1 | grep -oP 'v[0-9.]+')
     END_VERSIONS
     """
}

process filter_mutect {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}", pattern: "${meta.lib}*"

     input:
     tuple val(meta),
     path(vcf),
     path(tbi),
     path(stats),
     path(orientation_model)
     path genome
     path genome_fai
     path genome_dict

     output:
     tuple val(meta), path("${meta.lib}.Mutect2_filtered.vcf.gz"),     emit: mutect2_filtered_vcf
     tuple val(meta), path("${meta.lib}.Mutect2_filtered.vcf.gz.tbi"), emit: mutect2_filtered_vcf_tbi
     tuple val(meta), path("${meta.lib}.Mutect2_filtering.stats"),     emit: mutect2_filtering_stats
     path "versions.yml",                                               emit: versions

     stub:
     """
        touch "${meta.lib}.Mutect2_filtered.vcf.gz"
        touch "${meta.lib}.Mutect2_filtered.vcf.gz.tbi"
        touch "${meta.lib}.Mutect2_filtering.stats"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
     module load GATK/4.6.2.0
     gatk --java-options "-Xmx16g" FilterMutectCalls \
          -R ${genome} \
          -V ${vcf} \
          --stats ${stats} \
          --ob-priors ${orientation_model} \
          -O ${prefix}.Mutect2_filtered.vcf.gz \
          --filtering-stats ${prefix}.Mutect2_filtering.stats

     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
         gatk: \$(gatk --version 2>&1 | grep -oP 'v[0-9.]+')
     END_VERSIONS
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
