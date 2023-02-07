
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

//import workflows

include {Cutadapt} from './modules/cutadapt/cutadapt'
include {Fastqc} from './modules/qc/qc'
include {Multiqc} from './modules/qc/qc'
include {Mixcr_VCJtools} from './modules/misc/mixcr'
include {HLA_calls} from './workflows/hla_calls'
include {Star_rsem} from './workflows/star-rsem'
include {Fusion_calling} from './workflows/Fusion_calling'
include {Star_bam_processing} from './workflows/Star_bam_processing'
include {RNAseq_GATK} from './workflows/RNAseq_GATK'
include {QC_from_finalBAM} from './workflows/QC_from_finalBAM'
include {Annotation} from './workflows/Annotation'
include {QC_from_Star_bam} from './workflows/QC_from_Star_bam'

workflow {
    read_pairs              = Channel
                                .fromFilePairs(params.reads, flat: true)
                                .ifEmpty { exit 1, "Read pairs could not be found: ${params.reads}" }
    starfusion_db           = Channel.of(file(params.starfusion_db, checkIfExists:true))
    mixcr_license           = Channel.of(file(params.mixcr_license, checkIfExists:true))

// Trim away adapters
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

// QC with FastQC 
Fastqc(fqc_inputs.fqc_input)

if (params.run_upto_counts) {
 
  Star_rsem(Cutadapt.out)
}  else {

  Star_rsem(Cutadapt.out) 
  Starfusion_input = Star_rsem.out.star.flatMap{it -> [id: it[0], chimeric_junctions: it[4]]}
  Star_rsem.out.star.branch{ id, tbam, bam, bai, chimeric_junctions ->
              other: true
                  return( tuple(id, chimeric_junctions))} \
              .set{Starfusion_input_tmp}
  Starfusion_input = Starfusion_input_tmp.other.combine(starfusion_db)
  Fusion_calling (
         Cutadapt.out,
         Starfusion_input
     )
  PicardARG_input = Star_rsem.out.star.flatMap{it -> [id: it[0], chimeric_junctions: it[4]]}
  Star_rsem.out.star.branch{ id, tbam, bam, bai, chimeric_junctions ->
              other: true
                  return( tuple(id, bam, bai))} \
              .set{PicardARG_input}  
 // PicardARG_input.view()

  Star_bam_processing(PicardARG_input)
  HLA_calls(Cutadapt.out)
  QC_from_Star_bam(
      Star_bam_processing.out.picard_ARG,
      Star_bam_processing.out.picard_MD
  )
  RNAseq_GATK(Star_bam_processing.out.picard_MD)
  QC_from_finalBAM(RNAseq_GATK.out.GATK_RNAseq_bam)
  Annotation(
      RNAseq_GATK.out.SnpEff_vcf.combine(QC_from_finalBAM.out.hotspot_pileup, by:0),
      RNAseq_GATK.out.SnpEff_vcf
)
  }
  if (params.run_upto_counts) {

    multiqc_input = Fastqc.out.join(Star_rsem.out.star).join(Star_rsem.out.rsem)
    Multiqc(multiqc_input)
  }  else {  
    multiqc_input = Fastqc.out \
                       .join(QC_from_finalBAM.out.hotspot_pileup) \
                       .join(QC_from_finalBAM.out.coverageplot) \
                       .join(Star_rsem.out.star) \
                       .join(Star_rsem.out.rsem).join(QC_from_Star_bam.out.rnaseqc) \
                       .join(QC_from_Star_bam.out.circos) 
    multiqc_input.view() 
    Multiqc(multiqc_input)
  }

}
