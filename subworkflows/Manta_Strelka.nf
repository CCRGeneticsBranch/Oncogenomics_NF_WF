include {Manta} from '../modules/Variant_analysis/Manta_Strelka.nf'
include {Strelka} from '../modules/Variant_analysis/Manta_Strelka.nf'
include {Strelka_vcf_processing} from '../modules/Variant_analysis/Manta_Strelka.nf'
include {SnpEff} from '../modules/misc/snpEff'
include {Vcf2txt} from '../modules/misc/snpEff'

workflow Manta_Strelka {

    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
    strelka_config          = Channel.of(file(params.strelka_config, checkIfExists:true))
    dbNSFP2_4             = Channel.of(file(params.dbNSFP2_4, checkIfExists:true))
    dbNSFP2_4_tbi         = Channel.of(file(params.dbNSFP2_4_tbi, checkIfExists:true))
    Biowulf_snpEff_config  = Channel.of(file(params.Biowulf_snpEff_config, checkIfExists:true))
    strelka_snvch           = Channel.from("strelka.snvs")
take:
    tumor_bam
    normal_bam

main:
    Manta(
        tumor_bam,
        normal_bam.map { tuple -> tuple.take(tuple.size() - 1) },
        genome,
        genome_fai,
        genome_dict
    )


    Strelka(
        tumor_bam.combine(Manta.out.manta_indels_vcf,by:[0]).combine(Manta.out.manta_indels_tbi,by:[0]),
        normal_bam.map { tuple -> tuple.take(tuple.size() - 1) },
        genome,
        genome_fai,
        genome_dict,
        strelka_config
    )

    Strelka_vcf_processing(
        tumor_bam,
        Strelka.out.strelka_snv,
        Strelka.out.strelka_indels,
        normal_bam.map { tuple -> tuple.take(tuple.size() - 1) }
    )

    SnpEff(Strelka_vcf_processing.out.strelka_snv
               .combine(dbNSFP2_4)
               .combine(dbNSFP2_4_tbi)
               .combine(Biowulf_snpEff_config)
               .combine(strelka_snvch)
    )
    Vcf2txt(SnpEff.out.raw_snpeff.combine(strelka_snvch))


emit:
    strelka_indel_raw_vcf = Strelka_vcf_processing.out.strelka_indel
    strelka_snvs_raw_vcf = Strelka_vcf_processing.out.strelka_snv
    strelka_snpeff_snv_vcf2txt = Vcf2txt.out

}
