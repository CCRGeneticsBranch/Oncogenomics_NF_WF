
//include required modules
include {Picard_AddReadgroups} from '../modules/qc/picard'
include {Picard_MarkDuplicates} from '../modules/qc/picard'
include {Picard_CollectRNAseqmetrics} from '../modules/qc/picard'
include {Picard_CollectAlignmentSummaryMetrics} from '../modules/qc/picard'
include {RNAlibrary_customQC} from '../modules/qc/picard'

//Star_bam_processing workflow includes all the steps downstream of Star_RSEM workflow.
workflow Star_bam_processing {
     ref_flat                = Channel.of(file(params.ref_flat, checkIfExists:true))
     rRNA_interval           = Channel.of(file(params.rRNA_interval, checkIfExists:true))
     genome                  = Channel.of(file(params.genome, checkIfExists:true))
     genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
     genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))

    take: Coord_bam
          strandedness
          Fastqc_out
    main:
     Picard_AddReadgroups(Coord_bam)
     PicardCRS_input = Picard_AddReadgroups.out.combine(strandedness, by:[0]).combine(ref_flat).combine(rRNA_interval)
     Picard_CollectRNAseqmetrics(PicardCRS_input)

     Picard_CollectAlignmentSummaryMetrics(
         Picard_AddReadgroups.out
             .combine(genome)
     )

     Picard_MarkDuplicates(Picard_AddReadgroups.out)

     picard_markdup = Picard_MarkDuplicates.out.dedup_bam.join(Picard_MarkDuplicates.out.dedup_bam_bai)

     ch_versions = Picard_CollectRNAseqmetrics.out.versions.mix(Picard_MarkDuplicates.out.versions)

     RNAlib_qc_input = Picard_CollectRNAseqmetrics.out.rnaseq_metrics
                            .combine(Picard_CollectAlignmentSummaryMetrics.out, by:[0])
                            .combine(Fastqc_out, by:[0])
     RNAlibrary_customQC(RNAlib_qc_input)



    emit:

     picard_ARG = Picard_AddReadgroups.out
     picard_MD =  picard_markdup
     rnalib_custom_qc = RNAlibrary_customQC.out
     picard_rnaseqmetrics = Picard_CollectRNAseqmetrics.out.rnaseq_metrics
     picard_rnaseqmetrics_pdf = Picard_CollectRNAseqmetrics.out.rnaseq_metrics_pdf
     picard_alignmetrics = Picard_CollectAlignmentSummaryMetrics.out
     picard_version = ch_versions
}
