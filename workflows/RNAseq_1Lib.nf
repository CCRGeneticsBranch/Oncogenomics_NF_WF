include {Common_RNAseq_WF} from './Common_RNAseq_WF'
include {MakeHotSpotDB_1lib} from  '../modules/qc/plots'
include {Onelib_FormatInput} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {Multiqc} from '../modules/qc/qc'
include {Allstepscomplete} from '../modules/misc/Allstepscomplete'
include {AddAnnotation} from '../modules/annotation/annot'
workflow RNAseq_1lib {

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
    MakeHotSpotDB_1lib(Common_RNAseq_WF.out.pileup)
    format_input_ch = Common_RNAseq_WF.out.snpeff_vcf.combine(MakeHotSpotDB_1lib.out, by:[0])
    Onelib_FormatInput(format_input_ch)
    Annotation(
        Onelib_FormatInput.out)
    merged_ch = Common_RNAseq_WF.out.snpeff_vcf.combine(Annotation.out.rare_annotation,by:[0])
  
  AddAnnotation(merged_ch)   
    multiqc_input = Common_RNAseq_WF.out.Fastqc_out.join(Common_RNAseq_WF.out.pileup, by: [0])
                          .join(Common_RNAseq_WF.out.coverageplot, by: [0])
                          .join(Common_RNAseq_WF.out.chimeric_junction, by: [0])
                          .join(Common_RNAseq_WF.out.rsem_counts, by: [0]).join(Common_RNAseq_WF.out.rnaseqc, by: [0])
                          .join(Common_RNAseq_WF.out.circos_plot, by: [0])
    Multiqc(multiqc_input)
    final_inputs = Common_RNAseq_WF.out.fusion_calls.join(Annotation.out.final_annotation,by:[0]).join(Multiqc.out,by:[0])

    Allstepscomplete(final_inputs)

}