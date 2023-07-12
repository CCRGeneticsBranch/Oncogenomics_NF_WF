include {Exome_common_WF} from './Exome_common_WF.nf'
include {MakeHotSpotDB} from  '../modules/qc/plots'
include {FormatInput} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {AddAnnotation} from '../modules/annotation/annot'
include {CNVkitPooled} from '../modules/cnvkit/CNVkitPooled'
include {CNVkit_png} from '../modules/cnvkit/CNVkitPooled'

workflow Tumor_multiple_libs {

samples_exome = Channel.fromPath("Tumor_lib.csv")
.splitCsv(header:true)
.filter { row -> row.type == "Tumor" }
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

Exome_common_WF(samples_exome)

Exome_common_WF.out.pileup.map { meta, file ->
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

pileup_input_ch = combined_pileup_ch.map { tuple -> tuple.drop(1) } 
pileup_meta_ch = combined_pileup_ch.map { tuple -> tuple[0] }

MakeHotSpotDB(pileup_input_ch,
                   pileup_meta_ch
)

Exome_common_WF.out.HC_snpeff_snv_vcf2txt.map { meta, file ->
    meta2 = [
        id: meta.id,
        casename: meta.casename,
        type: meta.type,
        diagnosis: meta.diagnosis 
    ]
    [ meta2, file ]
  }.groupTuple()
   .map { meta, files -> [ meta, *files ] }
   .filter { tuple ->
    tuple.size() > 2
  }
   .set { combined_snpeff_ch }

formatinput_snpeff_ch  = combined_snpeff_ch.map { tuple -> tuple.drop(1) }

FormatInput(
    formatinput_snpeff_ch,
    MakeHotSpotDB.out
)
Annotation(FormatInput.out)

merged_ch = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.combine(Annotation.out.rare_annotation)
updated_tuples = merged_ch.map { tuple ->
    [tuple[0], tuple[1], tuple[3]]
}
AddAnnotation(updated_tuples)
cnvkit_clin_ex_v1 = Channel.of(file(params.cnvkit_clin_ex_v1, checkIfExists:true))

CNVkitPooled(
    Exome_common_WF.out.exome_final_bam.combine(cnvkit_clin_ex_v1)
)
CNVkit_png(CNVkitPooled.out.cnvkit_pdf)
}