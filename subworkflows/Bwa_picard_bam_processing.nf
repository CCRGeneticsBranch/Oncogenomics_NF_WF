include {BWA} from '../modules/mapping/bwa'
include {Picard_MarkDuplicates} from '../modules/qc/picard'
include {Fastqc} from '../modules/qc/qc'

workflow BWA_picard {

bwa_genomeindex   = Channel.of(file(params.bwa_genomeindex, checkIfExists:true))


take:samples_exome

main:

     //Initiate empty channel for versions
     ch_versions = Channel.empty()
     fastqc_input = samples_exome.map{ meta, r1, r2 -> [meta, [r1,r2]] }

     Fastqc(fastqc_input)

     ch_versions = ch_versions.mix(Fastqc.out.versions)

     BWA(samples_exome.combine(bwa_genomeindex))

     ch_versions = ch_versions.mix(BWA.out.versions)

     Picard_MarkDuplicates(BWA.out.bwa_bam.combine(BWA.out.bwa_bai,by:[0]))

     ch_versions = ch_versions.mix(Picard_MarkDuplicates.out.versions)


emit:
     bwa_bam = (BWA.out.bwa_bam.combine(BWA.out.bwa_bai,by:[0]))
     picard_MD =  Picard_MarkDuplicates.out.dedup_bam.combine(Picard_MarkDuplicates.out.dedup_bam_bai,by:[0])
     markdup_txt = Picard_MarkDuplicates.out.markdup_txt
     Fastqc_out = Fastqc.out.fastqc_results
     ch_versions = ch_versions

}
