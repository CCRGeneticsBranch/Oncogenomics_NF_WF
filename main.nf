
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
include {Arriba} from './modules/fusion/arriba'
include {Fusioncatcher} from './modules/fusion/fusioncatcher'
include {Starfusion} from './modules/fusion/starfusion'
include {Multiqc} from './modules/qc/qc'
include {Picard_AddReadgroups} from './modules/qc/picard'
include {Picard_MarkDuplicates} from './modules/qc/picard'
include {GATK_RNASeq_Trim} from './modules/RNAseq_GATK/GATK'
include {GATK_RNASeq_RTC_IR} from './modules/RNAseq_GATK/GATK'
include {GATK_RNASeq_BR_PR} from './modules/RNAseq_GATK/GATK'
include {Picard_CollectRNAseqmetrics} from './modules/qc/picard'
include {Picard_CollectAlignmentSummaryMetrics} from './modules/qc/picard'
include {Mixcr_VCJtools} from './modules/misc/mixcr'
include {HLAminer} from './modules/neoantigens/hlaminer'
include {Seq2HLA} from './modules/neoantigens/seq2hla.nf'
include {Hotspot_Coverage} from './modules/qc/plots.nf'
include {Hotspot_Boxplot} from './modules/qc/plots.nf'
include {Flagstat} from './modules/qc/plots.nf'
include {Bamutil} from './modules/qc/plots.nf'
include {RNAseq_HaplotypeCaller} from './modules/RNAseq_GATK/GATK'
include {MergeHLA} from './modules/neoantigens/mergeHLA.nf'
include {HotspotPileup} from './modules/qc/plots.nf'
include {SnpEff} from './modules/misc/snpEff'
include {Bam2tdf} from './modules/qc/plots.nf'
include {Vcf2txt} from './modules/misc/snpEff'
include {FormatInput} from './modules/annotation/annot'
include {Annovar} from './modules/annotation/annot'
//include {Genotyping} from  './modules/qc/qc'


workflow {
    read_pairs              = Channel
                                .fromFilePairs(params.reads, flat: true)
                                .ifEmpty { exit 1, "Read pairs could not be found: ${params.reads}" }

// Genome specifics
    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
    gtf                     = Channel.of(file(params.gtf, checkIfExists:true))
    chrom_sizes             = Channel.of(file(params.chrom_sizes, checkIfExists:true))

// STAR and RSEM
    star_genomeIndex        = Channel.of(file(params.star_genome_index, checkIfExists:true))
    rsemIndex               = Channel.of(file(params.rsem_index, checkIfExists:true))

// Arriba params
// These are now coming from the docker (ccbr_starplus)
    //blacklist               = Channel.of(file(params.blacklist), checkIfExists:true)
    //proteinDomains          = Channel.of(file(params.proteinDomains), checkIfExists:true)
    //cytobands               = Channel.of(file(params.cytobands), checkIfExists:true)

// Fusioncatcher db
    fusioncatcher_db        = Channel.of(file(params.fusioncatcher_db, checkIfExists:true))

// STARfusion db
    starfusion_db           = Channel.of(file(params.starfusion_db, checkIfExists:true))

// Mixcr license
    mixcr_license           = Channel.of(file(params.mixcr_license, checkIfExists:true))

// Picard and Genotyping
    // ref_flat                = Channel.of(file(params.ref_flat, checkIfExists:true))
     rRNA_interval           = Channel.of(file(params.rRNA_interval, checkIfExists:true))    
     phase1_1000g            = Channel.of(file(params.phase1_1000g, checkIfExists:true))
     Mills_and_1000g         = Channel.of(file(params.Mills_and_1000g, checkIfExists:true))
     Sites1000g4genotyping   = Channel.of(file(params.Sites1000g4genotyping, checkIfExists:true))
//     vcf2genotype            = Channel.of(file(params.vcf2genotype, checkIfExists:true))
//     vcf2loh                 = Channel.of(file(params.vcf2loh, checkIfExists:true))
     dbsnp                   = Channel.of(file(params.dbsnp, checkIfExists:true))
// hotspot bed files
     hg19_hotspot_pos         = Channel.of(file(params.hg19_hotspot_pos, checkIfExists:true))
     access_hotspot           = Channel.of(file(params.access_hotspot, checkIfExists:true))    
//  Annotation files
    annovar_data             = Channel.of(file(params.annovar_data, checkIfExists:true))    
    clinseq                  = Channel.of(file(params.clinseq, checkIfExists:true))
    cosmic                   = Channel.of(file(params.cosmic, checkIfExists:true))
    pcg                      = Channel.of(file(params.pcg, checkIfExists:true))
// biowulf snpEff config file
    Biowulf_snpEff_config = Channel.of(file(params.Biowulf_snpEff_config, checkIfExists:true))
    dbNSFP2_4             = Channel.of(file(params.dbNSFP2_4, checkIfExists:true))
    dbNSFP2_4_tbi         = Channel.of(file(params.dbNSFP2_4_tbi, checkIfExists:true))
// Trim away adapters
    Cutadapt(read_pairs)

    // combine raw fastqs and trimmed fastqs as input to fastqc
//    fastqc_input = Cutadapt.out.combine(read_pairs)
//    fastqc_input.branch { id1,trimr1,trimr2,id2,r1,r2 ->
//               fqc_input: id1 == id2
//                   return ( tuple (id1,r1,r2,trimr1,trimr2) )
//               other: true
//                   return ( tuple (id1,id2) )
//                  } \
//               .set { fqc_inputs }
    // fqc_inputs.fqc_input.view()

// QC with FastQC 
//    Fastqc(fqc_inputs.fqc_input)

// Align with STAR    
    Star(
        Cutadapt.out
            .combine(star_genomeIndex)
            .combine(gtf)
    )

// Count with RSEM
//    Rsem(
//        Star.out
//           .combine(rsemIndex)
//    )


// Fusion tools
// 1. Arriba
//    Arriba(
//        Cutadapt.out
//            .combine(genome)
//            .combine(star_genomeIndex)
//            .combine(gtf)
//    )
// 2. Fusioncatcher
//    Fusioncatcher(
//        Cutadapt.out
//            .combine(fusioncatcher_db)
//    )

// 3. Star-Fusion
 // Starfusion_input = Star.out.flatMap{it -> [id: it[0], chimeric_junctions: it[4]]}
//    Star.out.branch{ id, tbam, bam, bai, chimeric_junctions -> 
//                other: true
//                    return( tuple(id, chimeric_junctions))} \
//                .set{Starfusion_input_tmp}
//    Starfusion_input = Starfusion_input_tmp.other.combine(starfusion_db)                
    // Starfusion_input.view()
//    Starfusion(Starfusion_input)

//    Arriba.out
//        .combine(Fusioncatcher.out)
//        .combine(Starfusion.out)
//        .branch { id1, arriba_fusions_tsv, arriba_discarded_fusions_tsv, arriba_pdf, id2, fusioncatcher_final_list, fusioncatcher_summary, id3, starfusion_predictions_tsv ->
//            all_fusions: id1 == id2 == id3
//                return( tuple(id1, arriba_fusions_tsv, arriba_discarded_fusions_tsv, arriba_pdf, fusioncatcher_final_list, fusioncatcher_summary, starfusion_predictions_tsv))
//            other: true
//               return(tuple(id1,id2,id3))
//            } \
//            .set{merge_fusions_input}
//    merge_fusions_input.all_fusions.view()

// Mixcr
//    Mixcr_VCJtools (
//        Cutadapt.out
//           .combine(mixcr_license)
//    )

    // multiqc(fastqc.out)
     Picard_AddReadgroups(Star.out)    
    // Picard_CollectRNAseqmetrics(
    //Picard_AddReadgroups.out
    //         .combine(ref_flat)
    //         .combine(rRNA_interval) 
    // )    
    // Picard_CollectAlignmentSummaryMetrics(
    //     Picard_AddReadgroups.out
    //         .combine(genome)
    // )
     Picard_MarkDuplicates(Picard_AddReadgroups.out)


//    Genotyping(
//       Picard_AddReadgroups.out
//            .combine(genome)
//             .combine(genome_fai)
//             .combine(genome_dict)
//            .combine(Sites1000g4genotyping)
//    )

    GATK_RNASeq_Trim(
         Picard_MarkDuplicates.out
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
     )    
     GATK_RNASeq_RTC_IR(
         GATK_RNASeq_Trim.out
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
             .combine(phase1_1000g)
             .combine(Mills_and_1000g)
     )
     GATK_RNASeq_BR_PR(
         GATK_RNASeq_RTC_IR.out
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
             .combine(phase1_1000g)
             .combine(Mills_and_1000g)
     )

     Hotspot_Coverage(
        GATK_RNASeq_BR_PR.out
             .combine(chrom_sizes)
             .combine(access_hotspot)
     )

     Flagstat(GATK_RNASeq_BR_PR.out)

     Bamutil(
        GATK_RNASeq_BR_PR.out
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
     )

     HotspotPileup(
        GATK_RNASeq_BR_PR.out
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
             .combine(hg19_hotspot_pos)
     )

//     Bam2tdf(
//        GATK_RNASeq_BR_PR.out
//             .combine(genome)
//             .combine(genome_fai)
//             .combine(genome_dict)
//     )


//     Hotspot_Boxplot(Hotspot_Coverage.out)

// HLA prediction

  HLAminer(Cutadapt.out)

  Seq2HLA(Cutadapt.out)
  
  MergeHLA(
         Seq2HLA.out.combine(HLAminer.out, by:0)
  )

  RNAseq_HaplotypeCaller(
        GATK_RNASeq_BR_PR.out
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
             .combine(dbsnp)
     )

  SnpEff(
        RNAseq_HaplotypeCaller.out
             .combine(dbNSFP2_4)
             .combine(dbNSFP2_4_tbi)
     )
 Vcf2txt(SnpEff.out)

 FormatInput(Vcf2txt.out.combine(HotspotPileup.out, by:0))
 Annovar(
        FormatInput.out
             .combine(annovar_data)
             .combine(clinseq)
             .combine(pcg)
     )

}

