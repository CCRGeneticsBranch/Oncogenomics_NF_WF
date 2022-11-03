
// Using DSL-2
nextflow.enable.dsl=2


//Print out log
log.info """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         reads        : ${params.reads}
         """
         .stripIndent()

//import modules

include {cutadapt} from './modules/cutadapt/cutadapt'
include {fastqc} from './modules/qc/fastqc'
include {star} from './modules/mapping/star'
include {rsem} from './modules/quant/rsem'
include {multiqc} from './modules/qc/multiqc'
include {Picard_AddReadgroups} from './modules/qc/picard'
include {Picard_MarkDuplicates} from './modules/qc/picard'
include {GATK_RNASeq_Trim} from './modules/RNAseq_GATK/GATK'
include {GATK_RNASeq_RTC} from './modules/RNAseq_GATK/GATK'


// check input files

workflow{
    read_pairs = Channel
        .fromFilePairs(params.reads, flat: true)
        .ifEmpty { exit 1, "Read pairs could not be found: ${params.reads}" }
    star_genomeIndex = Channel.of(file(params.star_genome_index, checkIfExists:true))
    genome = Channel.of(file(params.genome, checkIfExists:true))
    genome_fai = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict = Channel.of(file(params.genome_dict, checkIfExists:true))
    gtf = Channel.of(file(params.gtf, checkIfExists:true))
    rsemIndex = Channel.of(file(params.rsem_index, checkIfExists:true))
    phase1_1000g = Channel.of(file(params.phase1_1000g, checkIfExists:true))
    Mills_and_1000g = Channel.of(file(params.Mills_and_1000g, checkIfExists:true))
    
    cutadapt(read_pairs)
    fastqc(cutadapt.out)
    star(
        cutadapt.out
            .combine(star_genomeIndex)
            .combine(gtf)
    )
    rsem(
        star.out
            .combine(rsemIndex)
    )
    multiqc(fastqc.out)
    Picard_AddReadgroups(star.out)    
    Picard_MarkDuplicates(Picard_AddReadgroups.out)
    GATK_RNASeq_Trim(
        Picard_MarkDuplicates.out
            .combine(genome)
            .combine(genome_fai)
            .combine(genome_dict)
    )    
    GATK_RNASeq_RTC(
        GATK_RNASeq_Trim.out
            .combine(genome)
            .combine(genome_fai)
            .combine(genome_dict)
            .combine(phase1_1000g)
            .combine(Mills_and_1000g)
    )

}

