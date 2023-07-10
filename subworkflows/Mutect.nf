include {Mutect} from '../modules/Variant_analysis/Mutect'
include {Mutect_order} from '../modules/Variant_analysis/Mutect'
include {SnpEff} from '../modules/misc/snpEff'
include {Vcf2txt} from '../modules/misc/snpEff'

workflow Mutect_WF {

    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
    dbsnp_138_b37_vcf       = Channel.of(file(params.dbsnp, checkIfExists:true))
    cosmic_v67_hg19_vcf     = Channel.of(file(params.cosmic_v67_hg19_vcf, checkIfExists:true))
    dbNSFP2_4             = Channel.of(file(params.dbNSFP2_4, checkIfExists:true))
    dbNSFP2_4_tbi         = Channel.of(file(params.dbNSFP2_4_tbi, checkIfExists:true))
    Biowulf_snpEff_config  = Channel.of(file(params.Biowulf_snpEff_config, checkIfExists:true))
    mutect_ch               = Channel.from("MuTect")
take:
    tumor_bam
    normal_bam

main:
     
    Mutect(
        tumor_bam,
        normal_bam.map { tuple -> tuple.take(tuple.size() - 1) },
        genome,
        genome_fai,
        genome_dict,
        dbsnp_138_b37_vcf,
        cosmic_v67_hg19_vcf
    )

    Mutect_order(Mutect.out.mutect_raw_vcf)

    SnpEff(Mutect_order.out
           .combine(dbNSFP2_4)
           .combine(dbNSFP2_4_tbi)
           .combine(Biowulf_snpEff_config)
           .combine(mutect_ch)
    )
    Vcf2txt(SnpEff.out.combine(mutect_ch))
emit:
    mutect_snpeff_snv_vcf2txt = Vcf2txt.out
    mutect_raw_vcf = Mutect_order.out
}