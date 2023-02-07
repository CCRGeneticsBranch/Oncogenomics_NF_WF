include {Picard_AddReadgroups} from '../modules/qc/picard'
include {Picard_MarkDuplicates} from '../modules/qc/picard'
include {Picard_CollectRNAseqmetrics} from '../modules/qc/picard'
include {Picard_CollectAlignmentSummaryMetrics} from '../modules/qc/picard'

workflow Star_bam_processing {
     ref_flat                = Channel.of(file(params.ref_flat, checkIfExists:true))
     rRNA_interval           = Channel.of(file(params.rRNA_interval, checkIfExists:true))
     genome                  = Channel.of(file(params.genome, checkIfExists:true))
     genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
     genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))

    take: Coord_bam
    main:
     Picard_AddReadgroups(Coord_bam)

     Picard_CollectRNAseqmetrics(
     Picard_AddReadgroups.out
             .combine(ref_flat)
             .combine(rRNA_interval)
     )
     Picard_CollectAlignmentSummaryMetrics(
         Picard_AddReadgroups.out
             .combine(genome)
     )
     Picard_MarkDuplicates(Picard_AddReadgroups.out)

    emit:
     
     picard_ARG = Picard_AddReadgroups.out
     picard_MD =  Picard_MarkDuplicates.out

}



