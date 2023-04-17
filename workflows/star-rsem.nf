
include {Star} from '../modules/mapping/star'
include {Rsem} from '../modules/quant/rsem'
include {Strandedness} from '../modules/quant/rsem'

workflow Star_rsem {

    star_genomeIndex        = Channel.of(file(params.star_genome_index, checkIfExists:true))
    rsemIndex               = Channel.of(file(params.rsem_index, checkIfExists:true))
    strandedness            = Channel.value(params.strandedness)
    gtf                     = Channel.of(file(params.gtf, checkIfExists:true))
    gtf_sorted               = Channel.of(file(params.gtf_sorted, checkIfExists:true))
    gtf_sorted_index         = Channel.of(file(params.gtf_sorted_index, checkIfExists:true))

    take: trimmed_fq

    main:
	Star(
            trimmed_fq
            .combine(star_genomeIndex)
            .combine(gtf)
    )
        Strandedness(
            Star.out
            .combine(gtf_sorted)
            .combine(gtf_sorted_index)
    )

        Rsem(Star.out.join(Strandedness.out, by: [0,1])
            .combine(rsemIndex)
    )
    emit:
     star =  Star.out
     rsem =  Rsem.out
     strandedness = Strandedness.out 
}

