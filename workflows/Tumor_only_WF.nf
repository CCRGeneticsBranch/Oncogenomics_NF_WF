include {Exome_common_WF} from './Exome_common_WF.nf'
include {MakeHotSpotDB} from  '../modules/qc/plots'
include {FormatInput} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {AddAnnotation} from '../modules/annotation/annot'
include {CNVkitPooled} from '../modules/cnvkit/CNVkitPooled'
include {CircosPlot} from '../modules/qc/qc'

workflow Tumor_only_WF {

samples_exome = Channel.fromPath(params.samplesheet1)
.splitCsv(header:true)
.filter { row -> row.type == "Tumor" }
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

Exome_common_WF(samples_exome)

pileup_input_ch = Exome_common_WF.out.pileup.map { tuple -> tuple[1] }    
pileup_meta_ch =Exome_common_WF.out.pileup.map { tuple -> tuple[0] }

MakeHotSpotDB(pileup_input_ch,
                       pileup_meta_ch
)

merged_loh_ch = Exome_common_WF.out.loh.map { tuple -> tuple[1] }
meta_merged_loh = Exome_common_WF.out.loh.map { tuple -> tuple[0] }

CircosPlot(
    merged_loh_ch,
    meta_merged_loh
)

formatinput_snpeff_ch = Exome_common_WF.out.snpeff_vcf.map { tuple -> tuple.drop(1) }

FormatInput(
        formatinput_snpeff_ch,
        MakeHotSpotDB.out)

Annotation(FormatInput.out)


merged_ch = Exome_common_WF.out.snpeff_vcf.combine(Annotation.out.rare_annotation,by:[0])
AddAnnotation(merged_ch)

cnvkit_clin_ex_v1 = Channel.of(file(params.cnvkit_clin_ex_v1, checkIfExists:true))

CNVkitPooled(
    Exome_common_WF.out.exome_final_bam.combine(cnvkit_clin_ex_v1)
)




}