include {BWA} from '../modules/mapping/bwa'
include {Picard_MarkDuplicates} from '../modules/qc/picard'

workflow BWA_picard {

bwa_genomeindex   = Channel.of(file(params.bwa_genomeindex, checkIfExists:true))


take:samples_exome

main:
     BWA(samples_exome.combine(bwa_genomeindex))
     Picard_MarkDuplicates(BWA.out)


emit: 
     bwa_bam = BWA.out
     picard_MD =  Picard_MarkDuplicates.out
}