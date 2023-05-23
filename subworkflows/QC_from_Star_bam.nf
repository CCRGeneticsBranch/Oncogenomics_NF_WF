include {RNAseQC} from '../modules/qc/qc'
include {CircosPlot} from  '../modules/qc/qc'
include {Genotyping} from  '../modules/qc/qc'

workflow QC_from_Star_bam {

   genome                  = Channel.of(file(params.genome, checkIfExists:true))
   genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
   genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
   gtf                     = Channel.of(file(params.gtf, checkIfExists:true))
   transcript_gtf          = Channel.of(file(params.transcript_gtf, checkIfExists:true))
   Sites1000g4genotyping   = Channel.of(file(params.Sites1000g4genotyping, checkIfExists:true))
   rRNA_interval           = Channel.of(file(params.rRNA_interval, checkIfExists:true))

   take:
         Picard_ARG_bam
         Picard_MD_bam
   main:
    Genotyping(
       Picard_ARG_bam
            .combine(Sites1000g4genotyping)
            .combine(genome)
            .combine(genome_fai)
            .combine(genome_dict)
    )

    CircosPlot(Genotyping.out)
     RNAseQC(
	Picard_MD_bam
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
             .combine(rRNA_interval)
             .combine(transcript_gtf)
     )
    emit:
         rnaseqc = RNAseQC.out
         circos = CircosPlot.out

}
