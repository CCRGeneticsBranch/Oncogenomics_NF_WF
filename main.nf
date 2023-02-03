
// Using DSL-2
nextflow.enable.dsl=2

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
include {Multiqc} from './modules/qc/qc'
include {Mixcr_VCJtools} from './modules/misc/mixcr'

include {FormatInput} from './modules/annotation/annot'
include {Annovar} from './modules/annotation/annot'
include {Custom_annotation} from './modules/annotation/annot'
include {Combine_annotation} from './modules/annotation/annot'

include {RNAseQC} from './modules/qc/qc'
include {CircosPlot} from  './modules/qc/qc'
include {Genotyping} from  './modules/qc/qc'

include {HLA_calls} from './workflows/hla_calls'
include {Star_rsem} from './workflows/star-rsem'
include {Fusion_calling} from './workflows/Fusion_calling'
include {Star_bam_processing} from './workflows/Star_bam_processing'
include {RNAseq_GATK} from './workflows/RNAseq_GATK'
include {QC_from_finalBAM} from './workflows/QC_from_finalBAM'

workflow {
    read_pairs              = Channel
                                .fromFilePairs(params.reads, flat: true)
                                .ifEmpty { exit 1, "Read pairs could not be found: ${params.reads}" }
/*
// Genome specifics
    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
    gtf                     = Channel.of(file(params.gtf, checkIfExists:true))
    transcript_gtf          = Channel.of(file(params.transcript_gtf, checkIfExists:true))
    chrom_sizes             = Channel.of(file(params.chrom_sizes, checkIfExists:true))


*/

// STARfusion db
    starfusion_db           = Channel.of(file(params.starfusion_db, checkIfExists:true))

// Mixcr license
//    mixcr_license           = Channel.of(file(params.mixcr_license, checkIfExists:true))

// Picard and Genotyping
     Sites1000g4genotyping   = Channel.of(file(params.Sites1000g4genotyping, checkIfExists:true))

//  Annotation files
    annovar_data             = Channel.of(file(params.annovar_data, checkIfExists:true))    
    clinseq                  = Channel.of(file(params.clinseq, checkIfExists:true))
    cosmic                   = Channel.of(file(params.cosmic, checkIfExists:true))
    pcg                      = Channel.of(file(params.pcg, checkIfExists:true))
    clinvar                  = Channel.of(file(params.clinvar, checkIfExists:true))
    hgmd                     = Channel.of(file(params.hgmd, checkIfExists:true))
    matchTrial               = Channel.of(file(params.matchTrial, checkIfExists:true))
    mcg                      = Channel.of(file(params.mcg, checkIfExists:true))
    DoCM                     = Channel.of(file(params.DoCM, checkIfExists:true))
    CanDL                    = Channel.of(file(params.CanDL, checkIfExists:true))
    targetted_cancer_care    = Channel.of(file(params.targetted_cancer_care, checkIfExists:true))
    civic                    = Channel.of(file(params.civic, checkIfExists:true))
    ACMG                     = Channel.of(file(params.ACMG, checkIfExists:true))
    hg19_BLsites             = Channel.of(file(params.hg19_BLsites, checkIfExists:true))           
    hg19_WLsites             = Channel.of(file(params.hg19_WLsites, checkIfExists:true))
// biowulf snpEff config file
    Biowulf_snpEff_config = Channel.of(file(params.Biowulf_snpEff_config, checkIfExists:true))
    dbNSFP2_4             = Channel.of(file(params.dbNSFP2_4, checkIfExists:true))
    dbNSFP2_4_tbi         = Channel.of(file(params.dbNSFP2_4_tbi, checkIfExists:true))

// Trim away adapters
    Cutadapt(read_pairs)
//    HLA_calls(Cutadapt.out)

    // combine raw fastqs and trimmed fastqs as input to fastqc
    fastqc_input = Cutadapt.out.combine(read_pairs)
    fastqc_input.branch { id1,trimr1,trimr2,id2,r1,r2 ->
               fqc_input: id1 == id2
                   return ( tuple (id1,r1,r2,trimr1,trimr2) )
               other: true
                   return ( tuple (id1,id2) )
                  } \
               .set { fqc_inputs }


// QC with FastQC 
Fastqc(fqc_inputs.fqc_input)

Star_rsem(Cutadapt.out)

/*
Run_upto_counts_only     = Channel.value('true')


if (Run_upto_counts_only == 'true') 
 
  workflow  {
  Cutadapt(read_pairs)
  fastqc_input = Cutadapt.out.combine(read_pairs)
  fastqc_input.branch { id1,trimr1,trimr2,id2,r1,r2 ->
             fqc_input: id1 == id2
                 return ( tuple (id1,r1,r2,trimr1,trimr2) )
             other: true
                 return ( tuple (id1,id2) )
                } \
             .set { fqc_inputs }

  Fastqc(fqc_inputs.fqc_input)
  Star_rsem(Cutadapt.out)
  }
else
  workflow {
  Cutadapt(read_pairs)
  Fastqc(fqc_inputs.fqc_input)
  Star_rsem(Cutadapt.out)
*/

  Starfusion_input = Star_rsem.out.Star.flatMap{it -> [id: it[0], chimeric_junctions: it[4]]}
  Star_rsem.out.Star.branch{ id, tbam, bam, bai, chimeric_junctions -> 
              other: true
                  return( tuple(id, chimeric_junctions))} \
              .set{Starfusion_input_tmp}
  Starfusion_input = Starfusion_input_tmp.other.combine(starfusion_db)                



  Fusion_calling (
         Cutadapt.out,
         Starfusion_input
     )

// Mixcr
//    Mixcr_VCJtools (
//        Cutadapt.out()
//           .combine(mixcr_license)
//    )

  PicardARG_input = Star_rsem.out.Star.flatMap{it -> [id: it[0], chimeric_junctions: it[4]]}
  Star_rsem.out.Star.branch{ id, tbam, bam, bai, chimeric_junctions ->
              other: true
                  return( tuple(id, bam, bai))} \
              .set{PicardARG_input}  
//  PicardARG_input.view()

  Star_bam_processing(PicardARG_input)
  RNAseq_GATK(Star_bam_processing.out.Picard_MD)
  QC_from_finalBAM(RNAseq_GATK.out.GATK_RNAseq_bam)

 


/*

    Genotyping(
       Picard_AddReadgroups.out
            .combine(Sites1000g4genotyping)
            .combine(genome)
            .combine(genome_fai)
            .combine(genome_dict)
    )

    CircosPlot(Genotyping.out)
     RNAseQC(
	Picard_MarkDuplicates.out
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
             .combine(rRNA_interval)
             .combine(transcript_gtf)
     )

//     multiqc(Hotspot_Boxplot.out)


 FormatInput(Vcf2txt.out.combine(HotspotPileup.out, by:0))
 Annovar(
        FormatInput.out
             .combine(annovar_data)
             .combine(clinseq)
             .combine(pcg)
     )
 Custom_annotation(
       FormatInput.out
             .combine(clinvar)
             .combine(hgmd)
             .combine(matchTrial)
             .combine(mcg)
             .combine(DoCM)
             .combine(CanDL)
             .combine(targetted_cancer_care)
             .combine(civic)
    )
 Combine_annotation(Annovar.out.combine(Custom_annotation.out, by:0)
             .combine(Vcf2txt.out, by:0)
             .combine(ACMG)
             .combine(hg19_BLsites)
             .combine(hg19_WLsites)
    )

*/
}


