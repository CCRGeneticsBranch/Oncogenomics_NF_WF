include {BWA} from '../modules/mapping/bwa'
include {Picard_MarkDuplicates} from '../modules/qc/picard'
include {Fastqc} from '../modules/qc/qc'

workflow BWA_picard {

bwa_genomeindex   = Channel.of(file(params.bwa_genomeindex, checkIfExists:true))


take:samples_exome

main:

     Fastqc(samples_exome.map { tuple -> tuple.drop(1) },
            samples_exome.map { tuple -> tuple[0] })
     BWA(samples_exome.combine(bwa_genomeindex))
     Picard_MarkDuplicates(BWA.out.bwa_bam.combine(BWA.out.bwa_bai,by:[0]))


emit:
     bwa_bam = (BWA.out.bwa_bam.combine(BWA.out.bwa_bai,by:[0]))
     picard_MD =  Picard_MarkDuplicates.out.dedup_bam.combine(Picard_MarkDuplicates.out.dedup_bam_bai,by:[0])
     markdup_txt = Picard_MarkDuplicates.out.markdup_txt
     Fastqc_out = Fastqc.out.fastqc_results

}
