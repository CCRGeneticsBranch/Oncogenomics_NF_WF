include {GATK_RTC_IR} from '../modules/RNAseq_GATK/GATK'
include {GATK_BR_PR} from '../modules/RNAseq_GATK/GATK'
include {Exome_HaplotypeCaller} from '../modules/RNAseq_GATK/GATK'
include {SnpEff} from '../modules/misc/snpEff'
include {Vcf2txt} from '../modules/misc/snpEff'

workflow Exome_GATK {
    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
    phase1_1000g            = Channel.of(file(params.phase1_1000g, checkIfExists:true))
    Mills_and_1000g         = Channel.of(file(params.Mills_and_1000g, checkIfExists:true))
    dbsnp                   = Channel.of(file(params.dbsnp, checkIfExists:true))
    dbNSFP2_4             = Channel.of(file(params.dbNSFP2_4, checkIfExists:true))
    dbNSFP2_4_tbi         = Channel.of(file(params.dbNSFP2_4_tbi, checkIfExists:true))
    Biowulf_snpEff_config  = Channel.of(file(params.Biowulf_snpEff_config, checkIfExists:true))
    HC_ch                  = Channel.from("HC")
take: MD_bam
      Capture_bed

main:
     GATK_RTC_IR(
       MD_bam
            .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
             .combine(phase1_1000g)
             .combine(Mills_and_1000g)
     )
     GATK_BR_PR(
        GATK_RTC_IR.out.Ir_bam.combine(GATK_RTC_IR.out.Ir_bai,by:[0])
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
             .combine(phase1_1000g)
             .combine(Mills_and_1000g)
     )
     Exome_HaplotypeCaller(
             GATK_BR_PR.out.final_bam.combine(GATK_BR_PR.out.final_bai,by:[0])
             .combine(Capture_bed,by:[0])
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
             .combine(dbsnp)
     )

     SnpEff(
        Exome_HaplotypeCaller.out.exome_HC
             .combine(dbNSFP2_4)
             .combine(dbNSFP2_4_tbi)
             .combine(Biowulf_snpEff_config)
             .combine(HC_ch)
     )
      Vcf2txt(SnpEff.out.raw_snpeff.combine(HC_ch))

emit:
     GATK_Exome_bam =  GATK_BR_PR.out.final_bam.combine(GATK_BR_PR.out.final_bai,by:[0])
     SnpEff_vcf      = SnpEff.out.raw_snpeff
     HC_snpeff_snv_vcf2txt = Vcf2txt.out
}
