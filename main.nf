// Using DSL-2
nextflow.enable.dsl=2

params.reads = "s3://ccr-genomics-testdata/testdata/Test*_R_T_R{1,2}.fastq.gz"


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
//include {star} from './modules/mapping/star'

workflow{
    read_pairs = Channel
        .fromFilePairs(params.reads, flat: true)
        .ifEmpty { exit 1, "Read pairs could not be found: ${params.reads}" }

    cutadapt(read_pairs)
    fastqc(cutadapt.out)
}

