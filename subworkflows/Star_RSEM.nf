
include {Star} from '../modules/mapping/star'
include {Rsem} from '../modules/quant/rsem'
include {Strandedness} from '../modules/quant/rsem'

workflow Star_RSEM {

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
            Star.out.genome_bam.join(Star.out.genome_bai, by: [0])
            .combine(gtf_sorted)
            .combine(gtf_sorted_index)
    )

        Rsem(Star.out.transcriptome_bam.join(Strandedness.out, by: [0])
            .combine(rsemIndex)
    )
    emit:
     transcriptome_bam =  Star.out.transcriptome_bam
     genome_bam = Star.out.genome_bam
     genome_bai = Star.out.genome_bai
     chimeric_junction = Star.out.chimeric_junction
     rsem_genes = Rsem.out.rsem_genes
     rsem_isoforms = Rsem.out.rsem_isoforms
     strandedness = Strandedness.out 
}

