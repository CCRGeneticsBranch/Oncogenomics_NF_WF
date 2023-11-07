include {Sequenza_utils} from '../modules/misc/Sequenza'
include {Sequenza} from '../modules/misc/Sequenza'
include {Sequenza_annot} from '../modules/misc/Sequenza'

workflow Sequenza_annotation {

    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
    gc50base                = Channel.of(file(params.gc50base, checkIfExists:true))
    sequenza_Rscript        = Channel.of(file(params.sequenza_Rscript, checkIfExists:true))
    combined_gene_list      = Channel.of(file(params.combined_gene_list, checkIfExists:true))

take:
    bam_variant_calling_pair
    tumor_target_capture

main:
    Sequenza_utils(
        bam_variant_calling_pair
            .combine(genome)
            .combine(genome_fai)
            .combine(genome_dict)
            .combine(gc50base)
    )
/*
    Sequenza(
        Sequenza_utils.out.sequenza_bin,
        sequenza_Rscript
    )


    Sequenza_annot(
    Sequenza.out.segments.combine(target_capture,by:[0]),
    combined_gene_list
    )

emit:
    sequenza = Sequenza_annot.out
    alternate = Sequenza.out.alternate
    versions = Sequenza_utils.out.versions
*/
}
