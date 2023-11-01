process Manta {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}",pattern: "candidateSmallIndels*"

     input:
     tuple val(meta),
     path(Nbam),
     path(Nindex),
     path(Tbam),
     path(Tindex),
     path(Tbed),
     path(genome),
     path(genome_fai),
     path(genome_dict)

     output:
     tuple val(meta),path("candidateSmallIndels.vcf.gz"), emit: manta_indels_vcf
     tuple val(meta),path("candidateSmallIndels.vcf.gz.tbi"), emit: manta_indels_tbi
     path "versions.yml"             , emit: versions

    stub:
    """
        touch "candidateSmallIndels.vcf.gz"

    """

    script:

    """
    configManta.py --normalBam ${Nbam} --tumorBam ${Tbam} --referenceFasta ${genome} --runDir ./manta
    ./manta/runWorkflow.py
    cp ./manta/results/variants/candidateSmallIndels.vcf.gz .
    cp ./manta/results/variants/candidateSmallIndels.vcf.gz.tbi .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        manta: \$(configManta.py 2>&1  |grep -E '^Version:'|sed 's/Version: //')
    END_VERSIONS
    """
}

process Strelka {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}",pattern: "somatic*"

     input:
     tuple val(meta),
     path(Nbam),
     path(Nindex),
     path(Tbam),
     path(Tindex),
     path(Tbed),
     path(vcf),
     path(tbi),
     path(genome),
     path(genome_fai),
     path(genome_dict),
     path(strelka_config)

     output:
     tuple val(meta),path("somatic.snvs.vcf.gz"), emit: strelka_snvs_raw_vcf
     tuple val(meta),path("somatic.snvs.vcf.gz.tbi"), emit: strelka_snvs_raw_tbi
     tuple val(meta),path("somatic.indels.vcf.gz"), emit: strelka_indels_vcf
     tuple val(meta),path("somatic.indels.vcf.gz.tbi"), emit: strelka_indels_tbi
     path "versions.yml"             , emit: versions
     stub:
     """
     touch "somatic.snvs.vcf.gz"
     touch "somatic.snvs.vcf.gz.tbi"
     touch "somatic.indels.vcf.gz"
     touch "somatic.indels.vcf.gz.tbi"

     """

    script:

    """
    configureStrelkaSomaticWorkflow.py --normalBam=${Nbam} --tumorBam=${Tbam} --referenceFasta=${genome} --config=${strelka_config} --indelCandidates ${vcf} --runDir=./strelka --exome
    ./strelka/runWorkflow.py -m local
    cp ./strelka/results/variants/somatic.snvs.vcf.gz .
    cp ./strelka/results/variants/somatic.snvs.vcf.gz.tbi .
    cp ./strelka/results/variants/somatic.indels.vcf.gz .
    cp ./strelka/results/variants/somatic.indels.vcf.gz.tbi .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        strelka: \$(configureStrelkaSomaticWorkflow.py 2>&1  |grep -E '^Version:'|sed 's/Version: //')
    END_VERSIONS
    """

}


process Strelka_vcf_processing {
     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

     input:
     tuple val(meta),
     path(Nbam),
     path(Nindex),
     path(Tbam),
     path(Tindex),
     path(Tbed),
     path(snvs_vcf),
     path(snvs_vcf_tbi),
     path(indels_vcf),
     path(indels_vcf_tbi)
    

     output:
     tuple val(meta),path("${meta.lib}.strelka.snvs.raw.vcf"), emit: strelka_snv
     tuple val(meta),path("${meta.lib}.strelka.indels.raw.vcf"), emit: strelka_indel
     path "versions.yml"             , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    vcftools --gzvcf ${snvs_vcf} --bed  ${Tbed} --out ${prefix}.strelka.snvs.raw.vcf --recode --keep-INFO-all
    mv -f ${prefix}.strelka.snvs.raw.vcf.recode.vcf ${prefix}.strelka.snvs.raw.vcf
    vcftools --gzvcf ${indels_vcf} --bed  ${Tbed} --out ${prefix}.strelka.indels.raw.vcf --recode --keep-INFO-all
    mv -f ${prefix}.strelka.indels.raw.vcf.recode.vcf ${prefix}.strelka.indels.raw.vcf
    sed -i "s/FORMAT\tNORMAL\tTUMOR/FORMAT\t${meta.normal_id}\t${prefix}/g" ${prefix}.strelka.snvs.raw.vcf
    sed -i "s/FORMAT\tNORMAL\tTUMOR/FORMAT\t${meta.normal_id}\t${prefix}/g" ${prefix}.strelka.indels.raw.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcftools: \$(vcftools 2>&1  |grep -E '^VCFtools'|sed 's/VCFtools //'|awk -F'[()]' '{print \$2}')
    END_VERSIONS
    """

}
