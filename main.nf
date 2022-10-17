// Using DSL-2
nextflow.enable.dsl=2

//params.reads = "s3://ccr-genomics-testdata/testdata/Test*_R_T_R{1,2}.fastq.gz"
params.genome_index = "s3://ccr-genomics-testdata/References/index-STAR_2.7.9a"
params.gtf = "s3://ccr-genomics-testdata/References/gencode.v37lift37.annotation_ERCC92.gtf"
params.rsem_index = "s3://ccr-genomics-testdata/References/rsem_1.3.2"
params.s3_bucket = "s3://agc-424336837382-us-east-1"
//reads_ch = Channel.fromFilePairs(params.reads)

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


workflow{
    read_pairs = Channel
        .fromFilePairs(params.reads, flat: true)
        .ifEmpty { exit 1, "Read pairs could not be found: ${params.reads}" }
    genomeIndex = Channel.of(file(params.genome_index, checkIfExists:true))
    gtf = Channel.of(file(params.gtf, checkIfExists:true))
    genomeIndex = Channel.of(file(params.genome_index, checkIfExists:true))
    rsemIndex = Channel.of(file(params.rsem_index, checkIfExists:true))
    cutadapt(read_pairs)
    fastqc(cutadapt.out)
    star(
        cutadapt.out
            .combine(genomeIndex)
            .combine(gtf)
    )
    rsem(
        star.out
            .combine(rsemIndex)
    )
    multiqc(fastqc.out)

}

