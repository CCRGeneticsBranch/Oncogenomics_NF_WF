include {Common_RNAseq_WF} from './Common_RNAseq_WF'
include {MakeHotSpotDB} from  '../modules/qc/plots'
include {FormatInput} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {Multiqc} from '../modules/qc/qc'
include {Allstepscomplete} from '../modules/misc/Allstepscomplete'
include {AddAnnotation} from '../modules/annotation/annot'
include {DBinput} from '../modules/misc/DBinput'
include {RNAqc_TrancriptCoverage} from '../modules/qc/picard'

workflow RNAseq_only {

samples_rnaseq = Channel.fromPath("RNAseq.csv")
.splitCsv(header:true)
.filter { row -> row.type == "RNAseq" }
.map { row ->
    def meta = [:]
    meta.id    =  row.sample
    meta.lib   =  row.library
    meta.sc    =  row.sample_captures
    meta.casename  = row.casename 
    meta.type     = row.type
    meta.diagnosis =row.Diagnosis
    def fastq_meta = []
    fastq_meta = [ meta,  file(row.read1), file(row.read2)  ]

    return fastq_meta
}


    Common_RNAseq_WF(samples_rnaseq)

rnalib_qc_list_ch = Common_RNAseq_WF.out.rnalib_custum_qc.map { tuple -> tuple[1] }
rnaseqmetrics_list_ch = Common_RNAseq_WF.out.picard_rnaseqmetrics.map { tuple -> tuple[1] }
rnaseqmetrics_meta_ch = Common_RNAseq_WF.out.picard_rnaseqmetrics.map { tuple -> tuple[0] }

RNAqc_TrancriptCoverage(
           rnalib_qc_list_ch,
           rnaseqmetrics_list_ch,
           rnaseqmetrics_meta_ch
)

    pileup_input_ch = Common_RNAseq_WF.out.pileup.map { tuple -> tuple[1] }
    pileup_meta_ch =Common_RNAseq_WF.out.pileup.map { tuple -> tuple[0] }

    MakeHotSpotDB(pileup_input_ch,
                  pileup_meta_ch
    )

    //MakeHotSpotDB(Common_RNAseq_WF.out.pileup)
  
    //format_input_ch = Common_RNAseq_WF.out.snpeff_vcf.combine(MakeHotSpotDB.out, by:[0])
    formatinput_snpeff_ch = Common_RNAseq_WF.out.snpeff_vcf.map { tuple -> tuple.drop(1) }


    FormatInput(
        formatinput_snpeff_ch,
        MakeHotSpotDB.out
    )
    Annotation(FormatInput.out)
    merged_ch = Common_RNAseq_WF.out.snpeff_vcf.combine(Annotation.out.rare_annotation,by:[0])
  
   AddAnnotation(merged_ch) 
   dbinput_snpeff_ch = Common_RNAseq_WF.out.snpeff_vcf.map{ tuple -> tuple.drop(1) }
   //dbinput_snpeff_ch.view()
   dbinput_annot_ch = AddAnnotation.out.map{ tuple -> tuple.drop(1) }
   //dbinput_annot_ch.view()
   dbinput_meta_ch = AddAnnotation.out.map { tuple -> tuple[0] }
   //dbinput_meta_ch.view()

   DBinput(
     dbinput_annot_ch,
     dbinput_snpeff_ch,
     dbinput_meta_ch
   )

    multiqc_input = Common_RNAseq_WF.out.Fastqc_out.join(Common_RNAseq_WF.out.pileup, by: [0])
                          .join(Common_RNAseq_WF.out.coverageplot, by: [0])
                          .join(Common_RNAseq_WF.out.chimeric_junction, by: [0])
                          .join(Common_RNAseq_WF.out.rsem_counts, by: [0]).join(Common_RNAseq_WF.out.rnaseqc, by: [0])
                          .join(Common_RNAseq_WF.out.circos_plot, by: [0])
    Multiqc(multiqc_input)
    final_inputs = Common_RNAseq_WF.out.fusion_calls.join(Annotation.out.final_annotation,by:[0]).join(Multiqc.out,by:[0])

    Allstepscomplete(final_inputs)

}