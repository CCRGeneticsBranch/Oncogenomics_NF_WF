include {Exome_common_WF} from './Exome_common_WF.nf'
include {MakeHotSpotDB} from  '../modules/qc/plots'
include {FormatInput} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {AddAnnotation} from '../modules/annotation/annot'
include {CNVkitPooled} from '../modules/cnvkit/CNVkitPooled'
include {CircosPlot} from '../modules/qc/qc'
include {DBinput} from '../modules/misc/DBinput'
include {Actionable_variants} from '../modules/Actionable'

workflow Exome_only_WF {

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

formatinput_input_ch = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.map{ meta, vcf -> [meta, [vcf]] }.join(MakeHotSpotDB.out)
FormatInput(formatinput_input_ch)

Annotation(FormatInput.out)

merged_ch = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.join(Annotation.out.rare_annotation,by:[0])
AddAnnotation(merged_ch)

dbinput_snpeff_ch = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.map{ meta, txt -> [meta, [txt]] }
dbinput_anno_ch = AddAnnotation.out.map{ meta, txt -> [meta, [txt]] }
dbinput_ch = dbinput_anno_ch.join(dbinput_snpeff_ch,by:[0])
DBinput(dbinput_ch)

cnvkit_input_bam = Exome_common_WF.out.exome_final_bam.branch{
        tumor: it[0].type == "tumor_DNA" }
        .map { tuple ->
        def meta = tuple[0]
        def bam = tuple[1]
        def bai = tuple[2]
        def cnv_ref = ''

        if (meta.sc == 'clin.ex.v1') {
            cnv_ref = params.cnvkit_clin_ex_v1
        } else if (meta.sc == 'clin.snv.v1') {
            cnv_ref = params.cnvkit_clin_snv_v1
        } else if (meta.sc == 'clin.cnv.v2') {
            cnv_ref = params.cnvkit_clin_cnv_v2
        } else if (meta.sc == 'clin.snv.v2') {
            cnv_ref = params.cnvkit_clin_snv_v2
        } else if (meta.sc == 'seqcapez.rms.v1') {
            cnv_ref = params.cnvkit_seqcapez_rms_v1
        } else if (meta.sc == 'seqcapez.hu.ex.v3') {
            cnv_ref = params.cnvkit_seqcapez_hu_ex_v3
        } else {
            return [meta, bam, bai]
        }
        return [meta, bam, bai, cnv_ref]
}
.filter { tuple ->
    tuple.size() == 4
}

cnvkit_input_bam.view()

/*


cnvkit_clin_ex_v1 = Channel.of(file(params.cnvkit_clin_ex_v1, checkIfExists:true))

CNVkitPooled(
    Exome_common_WF.out.exome_final_bam.combine(cnvkit_clin_ex_v1)
)



*/

}
