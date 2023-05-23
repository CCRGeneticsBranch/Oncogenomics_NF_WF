

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
include {test} from '../modules/test'


workflow RNASEQ {

samples_ch = Channel.fromPath(params.samplesheet)
.splitCsv(header:true)
.map { row ->
    def meta = [:]
    meta.id    =  row.sample
    meta.lib   =  row.library
    meta.sc    =  row.sample_captures
    def fastq_meta = []
    fastq_meta = [ meta,  file(row.read1), file(row.read2)  ]
    
    return fastq_meta
}
samples_ch.view()

test(samples_ch)
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

/*
  Starfusion_input = Star_rsem.out.star.flatMap{it -> [id: it[0], chimeric_junctions: it[5]]}
  Star_rsem.out.star.branch{ id, lib, tbam, bam, bai, chimeric_junctions ->
              other: true
                  return( tuple(id, lib, chimeric_junctions))} \
              .set{Starfusion_input_tmp}
     transcriptome_bam =  Star.out.transcriptome_bam
     genome_bam	= Star.out.genome_bam
     genome_bai	= Star.out.genome_bai
     chimeric_junction = Star.out.chimeric_junction
*/

  Starfusion_input = Star_rsem.out.chimeric_junction.combine(starfusion_db)

//  Starfusion_input = Starfusion_input_tmp.other.combine(starfusion_db)
  Fusion_calling (
         Cutadapt.out,
         Starfusion_input
     )

/*
  PicardARG_input = Star_rsem.out.star.flatMap{it -> [id: it[0], chimeric_junctions: it[5]]}
  Star_rsem.out.star.branch{ id, lib, tbam, bam, bai, chimeric_junctions ->
              other: true
                  return( tuple(id, lib, bam, bai))} \
              .set{PicardARG_input}
*/
  PicardARG_input = Star_rsem.out.genome_bam.combine(Star_rsem.out.genome_bai,by:[0,1])

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
                          .join(Star_rsem.out.chimeric_junction, by: [0, 1])
                          .join(Star_rsem.out.rsem, by: [0, 1]).join(QC_from_Star_bam.out.rnaseqc, by: [0, 1])
                          .join(QC_from_Star_bam.out.circos, by: [0, 1])
Multiqc(multiqc_input)

final_inputs = Fusion_calling.out.merge_fusion.join(Annotation.out.final_annotation,by:[0,1])

Allstepscomplete(final_inputs)

}

