include {GATK_RTC_IR} from '../modules/RNAseq_GATK/GATK'
include {GATK_BR_PR} from '../modules/RNAseq_GATK/GATK'

workflow Exome_GATK {
    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
    phase1_1000g            = Channel.of(file(params.phase1_1000g, checkIfExists:true))
    Mills_and_1000g         = Channel.of(file(params.Mills_and_1000g, checkIfExists:true))


take: MD_bam

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
        GATK_RTC_IR.out
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
             .combine(phase1_1000g)
             .combine(Mills_and_1000g)
     )
emit:
     GATK_Exome_bam =  GATK_BR_PR.out
}