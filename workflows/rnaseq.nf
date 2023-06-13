

//import modules
include {Cutadapt} from '../modules/cutadapt/cutadapt'
include {Fastqc} from '../modules/qc/qc'
include {Multiqc} from '../modules/qc/qc'
include {Mixcr_VCJtools} from '../modules/misc/mixcr'
include {Allstepscomplete} from '../modules/misc/Allstepscomplete'

// import subworkflows

include {HLA_calls} from '../subworkflows/hla_calls'
include {Star_rsem} from '../subworkflows/star-rsem'
include {Fusion_calling} from '../subworkflows/Fusion_calling'
include {Star_bam_processing} from '../subworkflows/Star_bam_processing'
include {RNAseq_GATK} from '../subworkflows/RNAseq_GATK'
include {QC_from_finalBAM} from '../subworkflows/QC_from_finalBAM'
include {Annotation} from '../subworkflows/Annotation'
include {QC_from_Star_bam} from '../subworkflows/QC_from_Star_bam'


workflow RNASEQ {
  
//params.sample_exome = "/data/khanlab/projects/Nextflow_dev/testing/exome_samplesheet.csv"
samples_rnaseq = Channel.fromPath(params.samplesheet1)
.splitCsv(header:true)
.filter { row -> row.type == "RNAseq" }
.map { row ->
    def meta = [:]
    meta.id    =  row.sample
    meta.lib   =  row.library
    meta.sc    =  row.sample_captures
    meta.casename  = row.casename 
    meta.type     = row.type
    def fastq_meta = []
    fastq_meta = [ meta,  file(row.read1), file(row.read2)  ]

    return fastq_meta
}
//.ifEmpty {exit 0, "No RNAseq samples found!" }
if (samples_rnaseq.count() == 0) {
    // Handle the case when samples_rnaseq channel is empty
    println("No RNAseq samples found!")
    
} else {
    // Proceed with the workflow using the non-empty samples_rnaseq channel

Cutadapt(samples_rnaseq)


//samples_rnaseq.view()


fastqc_input = Cutadapt.out.trim_reads.join(samples_rnaseq, by:[0])


  Fastqc(fastqc_input)

  Star_rsem(Cutadapt.out.trim_reads) 

  starfusion_db           = Channel.of(file(params.starfusion_db, checkIfExists:true))
  Starfusion_input = Star_rsem.out.chimeric_junction.combine(starfusion_db)

  Fusion_calling (
         Cutadapt.out,
         Starfusion_input
     )

  PicardARG_input = Star_rsem.out.genome_bam.combine(Star_rsem.out.genome_bai,by:[0])
  
  
  Star_bam_processing(
      PicardARG_input,
      Star_rsem.out.strandedness
  )



//  HLA_calls(Cutadapt.out)  docker needs to be fixed



  QC_from_Star_bam(
      Star_bam_processing.out.picard_ARG,
      Star_bam_processing.out.picard_MD
  )

  RNAseq_GATK(Star_bam_processing.out.picard_MD)

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

  Annotation(
      QC_from_finalBAM.out.hotspot_pileup,
      RNAseq_GATK.out.SnpEff_vcf
)

multiqc_input = Fastqc.out.join(QC_from_finalBAM.out.hotspot_pileup, by: [0])
                          .join(QC_from_finalBAM.out.coverageplot, by: [0])
                          .join(Star_rsem.out.chimeric_junction, by: [0])
                          .join(Star_rsem.out.rsem, by: [0]).join(QC_from_Star_bam.out.rnaseqc, by: [0])
                          .join(QC_from_Star_bam.out.circos, by: [0])
Multiqc(multiqc_input)

final_inputs = Fusion_calling.out.merge_fusion.join(Annotation.out.final_annotation,by:[0]).join(Multiqc.out,by:[0])

Allstepscomplete(final_inputs)

}

}

