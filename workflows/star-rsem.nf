
include {Star} from '../modules/mapping/star'
include {Rsem} from '../modules/quant/rsem'


workflow Star_rsem {

    star_genomeIndex        = Channel.of(file(params.star_genome_index, checkIfExists:true))
    rsemIndex               = Channel.of(file(params.rsem_index, checkIfExists:true))
    strandedness            = Channel.value(params.strandedness)
    gtf                     = Channel.of(file(params.gtf, checkIfExists:true))
    take: trimmed_fq
    main:
	Star(
            trimmed_fq
            .combine(star_genomeIndex)
            .combine(gtf)
    )
        Rsem(
	    Star.out
            .combine(rsemIndex)
            .combine(strandedness)
    )
    emit:
     star =  Star.out
     rsem =  Rsem.out
 
}

