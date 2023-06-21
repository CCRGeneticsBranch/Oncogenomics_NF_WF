include {Common_RNAseq_WF} from './Common_RNAseq_WF'
include {Lib2_RNAqc_TrancriptCoverage} from '../modules/qc/picard'
include {MakeHotSpotDB_2lib} from '../modules/qc/plots'
include {Twolib_FormatInput} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {AddAnnotation} from '../modules/annotation/annot'
workflow RNAseq_2lib {

samples_rnaseq = Channel.fromPath("RNA_lib.csv")
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
Common_RNAseq_WF.out.rnalib_custum_qc.map { meta, file ->
    meta2 = [
        id: meta.id,
        casename: meta.casename,
        type: meta.type
    ]
    [ meta2, file ]
  }.groupTuple()
   .map { meta, files -> [ meta, *files ] }
   .filter { tuple ->
    tuple.size() > 2
  }
   .set { combined_RNAqc_ch }

Common_RNAseq_WF.out.picard_rnaseqmetrics.map { meta, file ->
    meta2 = [
        id: meta.id,
        casename: meta.casename,
        type: meta.type
    ]
    [ meta2, file ]
  }.groupTuple()
   .map { meta, files -> [ meta, *files ] }
   .filter { tuple ->
    tuple.size() > 2
  }
   .set { combined_rnaseqmetrics_ch }


Lib2_RNAqc_TrancriptCoverage(
    combined_RNAqc_ch,
    combined_rnaseqmetrics_ch
)

Common_RNAseq_WF.out.pileup.map { meta, file ->
    meta2 = [
        id: meta.id,
        casename: meta.casename,
        type: meta.type
    ]
    [ meta2, file ]
  }.groupTuple()
   .map { meta, files -> [ meta, *files ] }
   .filter { tuple ->
    tuple.size() > 2
  }
   .set { combined_pileup_ch }

MakeHotSpotDB_2lib(combined_pileup_ch)
Common_RNAseq_WF.out.snpeff_vcf.map { meta, file ->
    meta2 = [
        id: meta.id,
        casename: meta.casename,
        type: meta.type
    ]
    [ meta2, file ]
  }.groupTuple()
   .map { meta, files -> [ meta, *files ] }
   .filter { tuple ->
    tuple.size() > 2
  }
   .set { combined_snpeff_ch }

Twolib_FormatInput(
    combined_snpeff_ch,
    MakeHotSpotDB_2lib.out
)
Annotation(Twolib_FormatInput.out)
//Common_RNAseq_WF.out.snpeff_vcf.view()

merged_ch = Common_RNAseq_WF.out.snpeff_vcf.combine(Annotation.out.rare_annotation)
updated_tuples = merged_ch.map { tuple ->
    [tuple[0], tuple[1], tuple[3]]
}
//updated_tuples.view()




AddAnnotation(updated_tuples)    

multiqc_input = Common_RNAseq_WF.out.Fastqc_out.join(Common_RNAseq_WF.out.pileup, by: [0])
                   .join(Common_RNAseq_WF.out.coverageplot, by: [0])
                   .join(Common_RNAseq_WF.out.chimeric_junction, by: [0])
                   .join(Common_RNAseq_WF.out.rsem_counts, by: [0]).join(Common_RNAseq_WF.out.rnaseqc, by: [0])
                   .join(Common_RNAseq_WF.out.circos_plot, by: [0])
//multiqc_input.view()

}
