
// Using DSL-2
nextflow.enable.dsl=2

// def randomstr = print new Random().with {(1..9).collect {(('a'..'z')).join()[ nextInt((('a'..'z')).join().length())]}.join()}

//Print out log
// log.info """\
//          R N A S E Q - N F   P I P E L I N E  
//          ===================================
//          git info     : $workflow.repository - $workflow.version [$workflow.commitId] --> gives error "Unexpected error [StackOverflowError]" ... reposity and commitId are null ... but version gives this error
//          manifest ver.: $workflow.manifest.version --> is null
//          """
log.info """\
         R N A S E Q - N F   P I P E L I N E  
         ===================================
         NF version   : $nextflow.version
         runName      : $workflow.runName
         username     : $workflow.userName
         configs      : $workflow.configFiles
         profile      : $workflow.profile
         cmd line     : $workflow.commandLine
         start time   : $workflow.start
         projectDir   : $workflow.projectDir
         launchDir    : $workflow.launchDir
         workdDir     : $workflow.workDir
         homeDir      : $workflow.homeDir
         reads        : ${params.reads}
         """
         .stripIndent()

//import modules

include {Cutadapt} from './modules/cutadapt/cutadapt'
include {Fastqc} from './modules/qc/qc'
include {Star} from './modules/mapping/star'
include {Rsem} from './modules/quant/rsem'
include {multiqc} from './modules/qc/qc'
include {Picard_AddReadgroups} from './modules/qc/picard'
include {Picard_MarkDuplicates} from './modules/qc/picard'
include {GATK_RNASeq_Trim} from './modules/RNAseq_GATK/GATK'
include {GATK_RNASeq_RTC_IR} from './modules/RNAseq_GATK/GATK'
include {GATK_RNASeq_BR_PR} from './modules/RNAseq_GATK/GATK'
include {Picard_CollectRNAseqmetrics} from './modules/qc/picard'
include {Picard_CollectAlignmentSummaryMetrics} from './modules/qc/picard'
// working on Genotyping process
// include {Genotyping} from  './modules/qc/qc'
// include {fusioncatcher} from './modules/fusion/fusion'


workflow {
    read_pairs              = Channel
                                .fromFilePairs(params.reads, flat: true)
                                .ifEmpty { exit 1, "Read pairs could not be found: ${params.reads}" }
    star_genomeIndex        = Channel.of(file(params.star_genome_index, checkIfExists:true))
    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
    gtf                     = Channel.of(file(params.gtf, checkIfExists:true))
    ref_flat                = Channel.of(file(params.ref_flat, checkIfExists:true))
    rRNA_interval           = Channel.of(file(params.rRNA_interval, checkIfExists:true))    
    rsemIndex               = Channel.of(file(params.rsem_index, checkIfExists:true))
    phase1_1000g            = Channel.of(file(params.phase1_1000g, checkIfExists:true))
    Mills_and_1000g         = Channel.of(file(params.Mills_and_1000g, checkIfExists:true))
    Sites1000g4genotyping   = Channel.of(file(params.Sites1000g4genotyping, checkIfExists:true))
    vcf2genotype            = Channel.of(file(params.vcf2genotype, checkIfExists:true))
    vcf2loh                 = Channel.of(file(params.vcf2loh, checkIfExists:true))
    
// Assigning inputs to all the process

    Cutadapt(read_pairs)

    // combine raw fastqs and trimmed fastqs as input to fastqc
    fastqc_input = Cutadapt.out.combine(read_pairs)
    fastqc_input.branch { id1,trimr1,trimr2,id2,r1,r2 ->
               fqc_input: id1 == id2
                   return ( tuple (id1,r1,r2,trimr1,trimr2) )
               other: true
                   return ( tuple (id1,id2) )
                  } \
               .set { fqc_inputs }
    fqc_inputs.fqc_input.view()
    
    Fastqc(fqc_inputs.fqc_input)
    
    Star(
       Cutadapt.out
           .combine(star_genomeIndex)
           .combine(gtf)
    )

    Rsem(
       Star.out
           .combine(rsemIndex)
    )

    // multiqc(fastqc.out)
    // Picard_AddReadgroups(star.out)    
    // Picard_CollectRNAseqmetrics(
    //     Picard_AddReadgroups.out
    //         .combine(ref_flat)
    //         .combine(rRNA_interval) 
    // )    
    // Picard_CollectAlignmentSummaryMetrics(
    //     Picard_AddReadgroups.out
    //         .combine(genome)
    // )
    // Picard_MarkDuplicates(Picard_AddReadgroups.out)
//    Genotyping(
//       Picard_AddReadgroups.out
//            .combine(genome)
//            .combine(Sites1000g4genotyping)
//            .combine(vcf2genotype)
//            .combine(vcf2loh)
//    )

    // GATK_RNASeq_Trim(
    //     Picard_MarkDuplicates.out
    //         .combine(genome)
    //         .combine(genome_fai)
    //         .combine(genome_dict)
    // )    
    // GATK_RNASeq_RTC_IR(
    //     GATK_RNASeq_Trim.out
    //         .combine(genome)
    //         .combine(genome_fai)
    //         .combine(genome_dict)
    //         .combine(phase1_1000g)
    //         .combine(Mills_and_1000g)
    // )
    // GATK_RNASeq_BR_PR(
    //     GATK_RNASeq_RTC_IR.out
    //         .combine(genome)
    //         .combine(genome_fai)
    //         .combine(genome_dict)
    //         .combine(phase1_1000g)
    //         .combine(Mills_and_1000g)
    // )

}

