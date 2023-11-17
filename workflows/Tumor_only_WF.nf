include {Exome_common_WF} from './Exome_common_WF.nf'
include {MakeHotSpotDB} from  '../modules/qc/plots'
include {FormatInput} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {AddAnnotation} from '../modules/annotation/annot'
include {CNVkitPooled} from '../modules/cnvkit/CNVkitPooled'
include {CircosPlot} from '../modules/qc/qc'
include {DBinput} from '../modules/misc/DBinput'
include {Actionable_variants} from '../modules/Actionable'

workflow Tumor_only_WF {

group               = Channel.from("variants")
   somatic_actionable_sites = Channel.of(file(params.somatic_actionable_sites, checkIfExists:true))
   combined_gene_list = Channel.of(file(params.combined_gene_list, checkIfExists:true))


samples_exome = Channel.fromPath("Exome.csv")
.splitCsv(header:true)
.filter { row -> row.type == "tumor_DNA" || row.type == "normal_DNA" }
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

MakeHotSpotDB_input = Exome_common_WF.out.pileup.map{ meta, pileup -> [meta, [pileup]] }

MakeHotSpotDB(MakeHotSpotDB_input)
Circosplot_input = Exome_common_WF.out.loh.map{ meta, loh -> [meta, [loh]] }
CircosPlot(Circosplot_input)
/*


formatinput_snpeff_ch = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.map { tuple -> tuple.drop(1) }

FormatInput(
        formatinput_snpeff_ch,
        MakeHotSpotDB.out)

Annotation(FormatInput.out)


merged_ch = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.combine(Annotation.out.rare_annotation,by:[0])
AddAnnotation(merged_ch)

dbinput_snpeff_ch = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.map{ tuple -> tuple.drop(1) }
dbinput_annot_ch = AddAnnotation.out.map{ tuple -> tuple.drop(1) }
dbinput_meta_ch = AddAnnotation.out.map { tuple -> tuple[0] }

DBinput(
    dbinput_annot_ch,
    dbinput_snpeff_ch,
    dbinput_meta_ch
)

Actionable_variants(DBinput.out
       .combine(Annotation.out.rare_annotation,by:[0])
       .combine(combined_gene_list)
       .combine(somatic_actionable_sites)
       .combine(group)
)


cnvkit_clin_ex_v1 = Channel.of(file(params.cnvkit_clin_ex_v1, checkIfExists:true))

CNVkitPooled(
    Exome_common_WF.out.exome_final_bam.combine(cnvkit_clin_ex_v1)
)


*/

}
