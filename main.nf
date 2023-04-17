
// Using DSL-2
nextflow.enable.dsl=2
import groovy.json.JsonSlurper

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
         casename     : ${params.casename}
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


// Reading the samplesheet

subject = Channel.fromPath(params.samplesheet)
          | splitCsv(header:true)
          | map { row-> tuple("${row.sample}","${row.library}","${row.read1}","${row.read2}", "${row.sample_captures}", "${row.Diagnosis}") }


//Creating tuple for read_pairs
pairs=    Channel.fromPath(params.samplesheet)
          | splitCsv(header:true)
          | map { row-> tuple("${row.sample}","${row.library}","${row.read1}","${row.read2}") }


//Creating target_captures_channel
capture_channel = subject.map { sample, library, read1, read2, sample_captures, diagnosis ->
    // set the hotspot file based on the capture type
    def target_file = ''
    if (sample_captures == 'access') {
        target_file = params.access_target
    } else if (sample_captures == 'polya_stranded') {
        target_file = params.polya_stranded_target
     } else if (sample_captures == 'polya') {
        target_file = params.polya_target
    }  else if (sample_captures == 'ribozero') {
        target_file = params.ribozero_target
    }  else if (sample_captures == 'SmartRNA') {
        target_file = params.smartrna_target
    }
    // return a tuple with the sample name and hotspot file
    tuple(sample, library, target_file)
}


// Trim away adapters
Cutadapt(pairs)


fastqc_input = Cutadapt.out.combine(pairs)
fastqc_input.branch { id1,library1,trimr1,trimr2,id2,library2,r1,r2 ->
           fqc_input: id1 == id2 & library1 == library2
               return ( tuple (id1,library1,r1,r2,trimr1,trimr2) )
           other: true
               return ( tuple (id1,id2) )
              } \
           .set { fqc_inputs }

// QC with FastQC 
Fastqc(fqc_inputs.fqc_input)

/*

if (params.run_upto_counts) {
 
  Star_rsem(Cutadapt.out)
}  else {

*/


  starfusion_db           = Channel.of(file(params.starfusion_db, checkIfExists:true))

  Star_rsem(Cutadapt.out) 


  Starfusion_input = Star_rsem.out.star.flatMap{it -> [id: it[0], chimeric_junctions: it[5]]}
  Star_rsem.out.star.branch{ id, lib, tbam, bam, bai, chimeric_junctions ->
              other: true
                  return( tuple(id, lib, chimeric_junctions))} \
              .set{Starfusion_input_tmp}
  Starfusion_input = Starfusion_input_tmp.other.combine(starfusion_db)
  Fusion_calling (
         Cutadapt.out,
         Starfusion_input
     )

  PicardARG_input = Star_rsem.out.star.flatMap{it -> [id: it[0], chimeric_junctions: it[5]]}
  Star_rsem.out.star.branch{ id, lib, tbam, bam, bai, chimeric_junctions ->
              other: true
                  return( tuple(id, lib, bam, bai))} \
              .set{PicardARG_input}

  Star_bam_processing(
      PicardARG_input,
      Star_rsem.out.strandedness
  )



//  HLA_calls(Cutadapt.out)



  QC_from_Star_bam(
      Star_bam_processing.out.picard_ARG,
      Star_bam_processing.out.picard_MD
  )



  RNAseq_GATK(Star_bam_processing.out.picard_MD)


  QC_from_finalBAM(
      RNAseq_GATK.out.GATK_RNAseq_bam,
      capture_channel
  )



  Annotation(
      RNAseq_GATK.out.SnpEff_vcf.combine(QC_from_finalBAM.out.hotspot_pileup, by:[0,1]),
      RNAseq_GATK.out.SnpEff_vcf
)

multiqc_input = Fastqc.out.join(QC_from_finalBAM.out.hotspot_pileup, by: [0, 1])
                          .join(QC_from_finalBAM.out.coverageplot, by: [0, 1])
                          .join(Star_rsem.out.star, by: [0, 1])
                          .join(Star_rsem.out.rsem, by: [0, 1]).join(QC_from_Star_bam.out.rnaseqc, by: [0, 1])
                          .join(QC_from_Star_bam.out.circos, by: [0, 1])
Multiqc(multiqc_input)

/*

  }
  if (params.run_upto_counts) {

    multiqc_input = Fastqc.out.join(Star_rsem.out.star).join(Star_rsem.out.rsem)
    Multiqc(multiqc_input)
  }  else {  

    multiqc_input = Fastqc.out \
                       .join(QC_from_finalBAM.out.hotspot_pileup) \
                       .join(QC_from_finalBAM.out.coverageplot) \
                       .join(Star_rsem.out.star) \
                       .join(Star_rsem.out.rsem).join(QC_from_Star_bam.out.rnaseqc) 
//                       .join(QC_from_Star_bam.out.circos) 
    Multiqc(multiqc_input)
  }

*/

}
