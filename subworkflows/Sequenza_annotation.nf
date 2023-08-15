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
    tumor_bam
    normal_bam 
    target_capture

main:
    Sequenza_utils(
        tumor_bam,
        normal_bam.map { tuple -> tuple.take(tuple.size() - 1) },
        genome,
        genome_fai,
        genome_dict,
        gc50base   
    )

    Sequenza(
        Sequenza_utils.out,
        sequenza_Rscript
    )
    
   
    Sequenza_annot(
    Sequenza.out.segments.combine(target_capture,by:[0]),
    combined_gene_list
    )

emit:
    sequenza = Sequenza_annot.out

}
