

//import modules
include {Cutadapt} from '../modules/cutadapt/cutadapt'
include {Kraken} from '../modules/qc/qc'
include {Krona} from '../modules/qc/qc'
include {Fastqc} from '../modules/qc/qc'
include {Fastq_screen} from '../modules/qc/qc'
include {Multiqc} from '../modules/qc/qc'
include {Mixcr_VCJtools} from '../modules/misc/mixcr'
include {Allstepscomplete} from '../modules/misc/Allstepscomplete'

// import subworkflows

include {HLA_calls} from '../subworkflows/hla_calls'
include {Star_RSEM} from '../subworkflows/Star_RSEM'
include {Fusion_calling} from '../subworkflows/Fusion_calling'
include {Star_bam_processing} from '../subworkflows/Star_bam_processing'
include {RNAseq_GATK} from '../subworkflows/RNAseq_GATK'
include {QC_from_finalBAM} from '../subworkflows/QC_from_finalBAM'
include {QC_from_Star_bam} from '../subworkflows/QC_from_Star_bam'



workflow Common_RNAseq_WF {

kraken_bacteria = Channel.of(file(params.kraken_bacteria, checkIfExists:true))
starfusion_db           = Channel.of(file(params.starfusion_db, checkIfExists:true))
fastq_screen_config         = Channel.of(file(params.fastq_screen_config, checkIfExists:true))

take:
     samples_rnaseq_ch

main:

ch_versions = Channel.empty()
Cutadapt(samples_rnaseq_ch)

ch_versions = ch_versions.mix(Cutadapt.out.versions)

Kraken(samples_rnaseq_ch
    .combine(kraken_bacteria)
)

ch_versions = ch_versions.mix(Kraken.out.versions)

Krona(Kraken.out.kraken_out)


fastqc_input = Cutadapt.out.trim_reads

  Fastqc(fastqc_input)

ch_versions = ch_versions.mix(Fastqc.out.versions)


  Fastq_screen_input = Cutadapt.out.trim_reads
                        .combine(fastq_screen_config)
                        .combine(Channel.fromPath(params.fastq_screen_db,type: 'dir', checkIfExists: true))

  Fastq_screen(Fastq_screen_input)



  Star_RSEM(Cutadapt.out.trim_reads)

 ch_versions = ch_versions.mix(Star_RSEM.out.star_version).mix(Star_RSEM.out.rsem_version)

  Starfusion_input = Star_RSEM.out.chimeric_junction.combine(starfusion_db)


  Fusion_calling (
         Cutadapt.out.trim_reads,
         Starfusion_input
     )

 ch_versions = ch_versions.mix(Fusion_calling.out.Arriba_version).mix(Fusion_calling.out.Fusioncatcher_version).mix(Fusion_calling.out.Starfusion_version)

  PicardARG_input = Star_RSEM.out.genome_bam.combine(Star_RSEM.out.genome_bai,by:[0])


  Star_bam_processing(
      PicardARG_input,
      Star_RSEM.out.strandedness,
      Fastqc.out.fastqc_results
  )

 ch_versions = ch_versions.mix(Star_bam_processing.out.picard_version)


HLA_calls(Cutadapt.out.trim_reads)

ch_version = ch_versions.mix(HLA_calls.out.version)

//Star_bam_processing.out.rnalib_custom_qc.view()
/*
Star_bam_processing.out.rnalib_custom_qc.map { meta, file ->
    new_meta = [
        id: meta.id,
        case: meta.case,
        type: meta.type
    ]
    [ new_meta, file ]
  }.groupTuple()
   .map { meta, files -> [ meta, *files ] }
   .filter { tuple ->
    tuple.size() > 2
  }
   .set { combinedChannel }


//combinedChannel.view()

*/

  QC_from_Star_bam(
      Star_bam_processing.out.picard_ARG,
      Star_bam_processing.out.picard_MD
  )

  RNAseq_GATK(Star_bam_processing.out.picard_MD)

  ch_versions = ch_versions.mix(RNAseq_GATK.out.GATK_version)

  capture_ch = RNAseq_GATK.out.GATK_RNAseq_bam
    .map { tuple ->
        def meta = tuple[0]
        def bam = tuple[1]
        def bai = tuple[2]
        def target_file = ''

        if (meta.sc == 'access') {
            target_file = params.access_target
        } else if (meta.sc == 'polya_stranded') {
            target_file = params.polya_stranded_target
        } else if (meta.sc == 'polya') {
            target_file = params.polya_target
        } else if (meta.sc == 'ribozero') {
            target_file = params.ribozero_target
        } else if (meta.sc == 'SmartRNA') {
            target_file = params.smartrna_target
        }

        return [meta,target_file]
    }



  QC_from_finalBAM(
      RNAseq_GATK.out.GATK_RNAseq_bam,
      capture_ch
  )

  ch_versions = ch_versions.mix(QC_from_finalBAM.out.flagstat_version).mix(QC_from_finalBAM.out.bamutil_version)

  emit:
  Fastqc_out = Fastqc.out.fastqc_results
  coverageplot = QC_from_finalBAM.out.coverageplot
  pileup = QC_from_finalBAM.out.hotspot_pileup
  snpeff_vcf = RNAseq_GATK.out.SnpEff_vcf
  chimeric_junction = Star_RSEM.out.chimeric_junction
  rsem_genes = Star_RSEM.out.rsem_genes
  rsem_isoforms = Star_RSEM.out.rsem_isoforms
  rnaseqc = QC_from_Star_bam.out.rnaseqc
  circos_plot = QC_from_Star_bam.out.circos
  fusion_calls = Fusion_calling.out.merge_fusion
  rnalib_custum_qc = Star_bam_processing.out.rnalib_custom_qc
  picard_rnaseqmetrics = Star_bam_processing.out.picard_rnaseqmetrics
  picard_rnaseqmetrics_pdf = Star_bam_processing.out.picard_rnaseqmetrics_pdf
  picard_alignmetrics = Star_bam_processing.out.picard_alignmetrics
  picard_MD = Star_bam_processing.out.picard_MD
  loh = QC_from_Star_bam.out.loh
  ch_versions = ch_versions
  gt = QC_from_Star_bam.out.gt
  strandedness = Star_RSEM.out.strandedness
  flagstat = QC_from_finalBAM.out.flagstat
  fastq_screen = Fastq_screen.out
  Bam = RNAseq_GATK.out.Bam
  hotspot_depth = QC_from_finalBAM.out.hotspot_depth
}
