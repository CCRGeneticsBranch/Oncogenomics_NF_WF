include {Exome_common_WF} from './Exome_common_WF.nf'
include {MakeHotSpotDB
        CoveragePlot} from  '../modules/qc/plots'
include {FormatInput} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {AddAnnotation} from '../modules/annotation/annot'
include {CNVkitPooled
        CNVkit_png} from '../modules/cnvkit/CNVkitPooled'
include {CircosPlot
        Genotyping_Sample
        Multiqc
        QC_summary_Patientlevel} from '../modules/qc/qc'
include {DBinput} from '../modules/misc/DBinput'
include {Actionable_variants} from '../modules/Actionable'
include {TcellExtrect} from '../modules/misc/TcellExtrect'
include {CUSTOM_DUMPSOFTWAREVERSIONS} from '../modules/nf-core/dumpsoftwareversions/main.nf'

workflow Exome_only_WF {

   somatic_actionable_sites = Channel.of(file(params.somatic_actionable_sites, checkIfExists:true))
   combined_gene_list = Channel.of(file(params.combined_gene_list, checkIfExists:true))
   genome_version_tcellextrect         = Channel.of(params.genome_version_tcellextrect)
   Pipeline_version = Channel.from(params.Pipeline_version)


samples_exome = Channel.fromPath("Exome.csv")
.splitCsv(header:true)
.filter { row -> row.type == "tumor_DNA" || row.type == "normal_DNA" || row.type == "cell_line_DNA" }
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


genotyping_input = (Exome_common_WF.out.gt.map{ meta, gt -> [meta, [gt]] })
Genotyping_Sample(genotyping_input,
                Pipeline_version)
Combined_coverage = Exome_common_WF.out.coverage.map{meta, coverage -> [meta, [coverage] ] }
CoveragePlot(Combined_coverage)

formatinput_input_ch = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.map{ meta, vcf -> [meta, [vcf]] }.join(MakeHotSpotDB.out)
FormatInput(formatinput_input_ch)

Annotation(FormatInput.out)

ch_versions = Exome_common_WF.out.ch_versions.mix(Annotation.out.version)

merged_ch = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.join(Annotation.out.rare_annotation,by:[0])
AddAnnotation(merged_ch)

dbinput_snpeff_ch = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.map{ meta, txt -> [meta, [txt]] }
dbinput_anno_ch = AddAnnotation.out.map{ meta, txt -> [meta, [txt]] }
dbinput_ch = dbinput_anno_ch.join(dbinput_snpeff_ch,by:[0])
DBinput(dbinput_ch)

cnvkit_input_bam = Exome_common_WF.out.exome_final_bam.branch{
        tumor: it[0].type == "tumor_DNA" || it[0].type == "cell_line_DNA"  }
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
        } else if (meta.sc == 'agilent.v7') {
            cnv_ref = params.cnvkit_agilent.v7
        } else {
            return [meta, bam, bai]
        }
        return [meta, bam, bai, cnv_ref]
}
.filter { tuple ->
    tuple.size() == 4
}

//cnvkit_input_bam.view()
cnvkit_input_bam|CNVkitPooled

CNVkitPooled.out.cnvkit_pdf|CNVkit_png

TcellExtrect(
    Exome_common_WF.out.exome_final_bam
    .join(Exome_common_WF.out.target_capture_ch,by:[0])
    .combine(genome_version_tcellextrect)
)

QC_summary_Patientlevel(Exome_common_WF.out.exome_qc)

multiqc_input = Exome_common_WF.out.Fastqc_out
    .join(Exome_common_WF.out.pileup, by: [0])
    .join(Exome_common_WF.out.kraken,by: [0])
    .join(Exome_common_WF.out.verifybamid,by: [0])
    .join(Exome_common_WF.out.hsmetrics,by: [0])
    .join(TcellExtrect.out,by: [0])
    .join(Exome_common_WF.out.fastq_screen,by: [0])
    .join(Exome_common_WF.out.flagstat,by: [0])
    .join(Exome_common_WF.out.markdup_txt,by: [0])
    .join(Exome_common_WF.out.krona,by: [0])

multiqc_input_ch = multiqc_input.map{ files ->
    if (files instanceof List) {
        return [files[0], files[1..-1]]
    } else {
        return files
    }
}

Multiqc(multiqc_input_ch)

ch_versions = ch_versions.mix(Multiqc.out.versions)
combine_versions  = ch_versions.unique().collectFile(name: 'collated_versions.yml')
custom_versions_input = Multiqc.out.multiqc_report
        .combine(combine_versions).map{ meta, multiqc, version -> [meta, version] }
        .combine(Pipeline_version)

CUSTOM_DUMPSOFTWAREVERSIONS(custom_versions_input)

}
