process Manta {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(Tbam),path(Tindex),path(Tbed)
     tuple val(meta2),path(Nbam),path(Nindex)
     path genome
     path genome_fai
     path genome_dict
     
     output:
     tuple val(meta),
     path("candidateSmallIndels.vcf.gz"),
     path("candidateSmallIndels.vcf.gz.tbi")

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
    """
}

process Strelka {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(Tbam),path(Tindex),path(Tbed),path(vcf),path(tbi)
     tuple val(meta2),path(Nbam),path(Nindex)
     path genome
     path genome_fai
     path genome_dict
     path strelka_config

     output:
     tuple val(meta),
     path("somatic.snvs.vcf.gz"),
     path("somatic.snvs.vcf.gz.tbi"),
     path("somatic.indels.vcf.gz"),
     path("somatic.indels.vcf.gz.tbi")

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
    
    """

}


process Strelka_vcf_processing {
     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),path(Tbam),path(Tindex),path(Tbed)
     tuple val(meta),path(snvs_vcf),path(snvs_vcf_tbi),path(indels_vcf),path(indels_vcf_tbi)
     tuple val(meta2),path(Nbam),path(Nindex)

     output:
     tuple val(meta),
     path("${meta.lib}.strelka.snvs.raw.vcf"),
     path("${meta.lib}.strelka.indels.raw.vcf")
    
    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    vcftools --gzvcf ${snvs_vcf} --bed  ${Tbed} --out ${prefix}.strelka.snvs.raw.vcf --recode --keep-INFO-all
    mv -f ${prefix}.strelka.snvs.raw.vcf.recode.vcf ${prefix}.strelka.snvs.raw.vcf
    vcftools --gzvcf ${indels_vcf} --bed  ${Tbed} --out ${prefix}.strelka.indels.raw.vcf --recode --keep-INFO-all
    mv -f ${prefix}.strelka.indels.raw.vcf.recode.vcf ${prefix}.strelka.indels.raw.vcf
    sed -i "s/FORMAT\tNORMAL\tTUMOR/FORMAT\t${meta2.lib}\t${prefix}/g" ${prefix}.strelka.snvs.raw.vcf
    sed -i "s/FORMAT\tNORMAL\tTUMOR/FORMAT\t${meta2.lib}\t${prefix}/g" ${prefix}.strelka.indels.raw.vcf

    """

}