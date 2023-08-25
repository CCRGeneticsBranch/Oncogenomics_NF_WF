
include {Arriba} from '../modules/fusion/arriba'
include {Fusioncatcher} from '../modules/fusion/fusioncatcher'
include {Starfusion} from '../modules/fusion/starfusion'
include {Mergefusion} from '../modules/fusion/merge'
workflow Fusion_calling {
    fusioncatcher_db        = Channel.of(file(params.fusioncatcher_db, checkIfExists:true))
    starfusion_db           = Channel.of(file(params.starfusion_db, checkIfExists:true))
    star_genomeIndex        = Channel.of(file(params.star_genome_index, checkIfExists:true))
    rsemIndex               = Channel.of(file(params.rsem_index, checkIfExists:true))
    strandedness            = Channel.value(params.strandedness)
    gtf                     = Channel.of(file(params.gtf, checkIfExists:true))
    genome                  = Channel.of(file(params.genome, checkIfExists:true))

    take:
          trimmed_fq
          starfusioninput
    main:
    Arriba(
	trimmed_fq
            .combine(genome)
            .combine(star_genomeIndex)
            .combine(gtf)
    )
    Fusioncatcher(
        trimmed_fq
            .combine(fusioncatcher_db)
    )
    Starfusion(starfusioninput)

    Mfinput = Arriba.out.arriba_fusion.join(Fusioncatcher.out.fc_output, by: [0]).join(Starfusion.out,by: [0])
    Mergefusion(Mfinput)
  

    emit:
     merge_fusion =  Mergefusion.out
 
}
