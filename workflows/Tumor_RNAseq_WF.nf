include {Common_RNAseq_WF} from './Common_RNAseq_WF'
include {Exome_common_WF} from './Exome_common_WF.nf'
include {MakeHotSpotDB
        Hotspot_Boxplot} from '../modules/qc/plots'
include {Exome_QC} from '../modules/qc/qc.nf'
include {Vcf2txt} from '../modules/misc/snpEff'
include {FormatInput
        AddAnnotation_TN} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {DBinput} from '../modules/misc/DBinput'
include {CNVkitPooled
        CNVkit_png} from '../modules/cnvkit/CNVkitPooled'
include {QC_summary_Patientlevel} from '../modules/qc/qc'
include {Genotyping_Sample
        Multiqc
        CircosPlot} from '../modules/qc/qc'
include {Actionable_fusion} from '../modules/Actionable.nf'
include {Fusion_Annotation
        Merge_fusion_annotation} from '../modules/annotation/Fusion_Annotation'
include {Combine_customRNAQC
        RNAqc_TrancriptCoverage} from '../modules/qc/picard'
include {CUSTOM_DUMPSOFTWAREVERSIONS} from '../modules/nf-core/dumpsoftwareversions/main.nf'


def combine_exome_rnaseq_libraries = { exomelib, RNAlib ->
    exomelib.cross(RNAlib).map { exome, rnaseq ->
                def meta = exome[1]
                [
                    meta + [
                        rna_lib: rnaseq[1].lib,
                        rna_type: rnaseq[1].type,
                        RNA_sc: rnaseq[1].sc
                    ],

                    [exome[2], rnaseq[2]]
                ]
            }
}



workflow Tumor_RNAseq_WF {

    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
    pfamdb  = Channel.of(file(params.pfamdb, checkIfExists:true))
    genome_version_fusion_annotation =  Channel.from(params.genome_version_fusion_annotation)
    genome_version = Channel.from(params.genome_version)
    Pipeline_version = Channel.from(params.Pipeline_version)


// Parse the samplesheet to generate fastq tuples
samples = Channel.fromPath("Tumor_RNAseq.csv")
.splitCsv(header:true)
.filter { row -> row.type == "tumor_DNA" || row.type == "cell_line_DNA" || row.type == "tumor_RNA" }
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

samples_branch = samples.branch{
        exome: it[0].type == "tumor_DNA" || it[0].type == "cell_line_DNA"
        rnaseq: it[0].type == "tumor_RNA"
}

samples_branch.rnaseq|Common_RNAseq_WF
samples_branch.exome|Exome_common_WF

ch_versions = Exome_common_WF.out.ch_versions.mix(Common_RNAseq_WF.out.ch_versions)



//RNA library pileup in  [meta.id, meta, file] format
pileup_samples_rnaseq_to_cross = Common_RNAseq_WF.out.pileup.map{ meta, pileup -> [ meta.id, meta, pileup ] }

//Tumor library pileup in [meta.id, meta, file] format
pileup_samples_tumor_to_cross = Exome_common_WF.out.pileup.map{ meta, pileup -> [ meta.id, meta, pileup ] }

pileup_pair = combine_exome_rnaseq_libraries(pileup_samples_tumor_to_cross,pileup_samples_rnaseq_to_cross)

MakeHotSpotDB(pileup_pair)

//Tumor sample vcftxt  in [meta.id, meta, file] format
HC_snpeff_snv_vcftxt_samples_tumor_to_cross = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.map{ meta, snpeff_snv_vcftxt -> [ meta.id, meta, snpeff_snv_vcftxt ] }

//RNAseq sample vcftxt  in [meta.id, meta, file] format
HC_snpeff_snv_vcftxt_samples_rna_to_cross = Common_RNAseq_WF.out.snpeff_vcf.map{ meta, snpeff_snv_vcftxt -> [ meta.id, meta, snpeff_snv_vcftxt ] }

Combined_snpeff_vcf2txt_ch = combine_exome_rnaseq_libraries(HC_snpeff_snv_vcftxt_samples_tumor_to_cross,HC_snpeff_snv_vcftxt_samples_rna_to_cross)

format_input_ch = Combined_snpeff_vcf2txt_ch.join(MakeHotSpotDB.out,by:[0])

FormatInput(format_input_ch)

Annotation(FormatInput.out)

ch_versions = ch_versions.mix(Annotation.out.version)

addannotation_input_ch = Combined_snpeff_vcf2txt_ch.join(Annotation.out.rare_annotation,by:[0])

AddAnnotation_TN(addannotation_input_ch)

Combined_AddAnnotation_TN_ch = AddAnnotation_TN.out.Tumor_hc_anno_txt
            .join(AddAnnotation_TN.out.RNA_hc_anno_txt,by:[0])
            .map{meta, tumor, rnaseq -> [meta, [tumor, rnaseq]]}
Combined_dbinput_ch = Combined_AddAnnotation_TN_ch.join(Combined_snpeff_vcf2txt_ch,by:[0])
DBinput(Combined_dbinput_ch)

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

cnvkit_input_bam|CNVkitPooled

CNVkitPooled.out.cnvkit_pdf|CNVkit_png

exome_hotspot_depth_status_tumor_to_cross = Exome_common_WF.out.hotspot_depth.map{meta, tumor -> [ meta.id, meta, tumor ] }
rnaseq_hotspot_depth = Common_RNAseq_WF.out.hotspot_depth.map{ meta, rnaseq -> [ meta.id, meta, rnaseq ] }
Combined_hotspot_ch = combine_exome_rnaseq_libraries(exome_hotspot_depth_status_tumor_to_cross,rnaseq_hotspot_depth)
Hotspot_Boxplot(Combined_hotspot_ch)

exome_genotyping_status_tumor_to_cross = Exome_common_WF.out.gt.map{ meta, tumor -> [ meta.id, meta, tumor ] }
rnaseq_genotyping_to_cross = Common_RNAseq_WF.out.gt.map{ meta, gt -> [ meta.id, meta, gt ] }
Combined_genotyping_ch = combine_exome_rnaseq_libraries(exome_genotyping_status_tumor_to_cross,rnaseq_genotyping_to_cross)
Genotyping_Sample(Combined_genotyping_ch)

exome_loh_status_tumor_to_cross = Exome_common_WF.out.loh.map{ meta, tumor -> [ meta.id, meta, tumor ] }
rnaseq_loh_to_cross = Common_RNAseq_WF.out.loh.map{ meta, loh -> [ meta.id, meta, loh ] }
Combined_loh_ch = combine_exome_rnaseq_libraries(exome_loh_status_tumor_to_cross,rnaseq_loh_to_cross)
CircosPlot(Combined_loh_ch)


 //RNA lib processing steps
actionable_fusion_input = Common_RNAseq_WF.out.fusion_calls.map{ meta, fusion -> [meta, [fusion]] }
Actionable_fusion(actionable_fusion_input)

Fusion_Annotation_input = Common_RNAseq_WF.out.rsem_isoforms
                        .join(Common_RNAseq_WF.out.fusion_calls, by:[0])
                        .combine(pfamdb)
                        .combine(genome)
                        .combine(genome_version_fusion_annotation)
                        .combine(genome_version)
Fusion_Annotation(Fusion_Annotation_input)

merge_fusion_anno_input = Fusion_Annotation.out.map{ meta, fusion -> [meta, [fusion]] }

Merge_fusion_annotation(merge_fusion_anno_input.combine(genome_version))

Combine_customRNAQC(Common_RNAseq_WF.out.rnalib_custum_qc.map{ meta, qc -> [meta, [qc]] })

RNAqc_TrancriptCoverage(Common_RNAseq_WF.out.picard_rnaseqmetrics.map{ meta, qc -> [meta, [qc]] })

//DNA lib processing

QC_summary_Patientlevel(Exome_common_WF.out.exome_qc)

multiqc_rnaseq_input = Common_RNAseq_WF.out.Fastqc_out.join(Common_RNAseq_WF.out.pileup, by: [0])
                      .join(Common_RNAseq_WF.out.coverageplot, by: [0])
                      .join(Common_RNAseq_WF.out.chimeric_junction, by: [0])
                      .join(Common_RNAseq_WF.out.rsem_genes, by: [0])
                      .join(Common_RNAseq_WF.out.rnaseqc, by: [0])
                      .join(Common_RNAseq_WF.out.circos_plot, by: [0])
                      .join(Common_RNAseq_WF.out.strandedness, by: [0])
                      .join(Common_RNAseq_WF.out.rnalib_custum_qc, by: [0])
                      .join(Common_RNAseq_WF.out.picard_rnaseqmetrics, by: [0])
                      .join(Common_RNAseq_WF.out.picard_rnaseqmetrics_pdf, by: [0])
                      .join(Common_RNAseq_WF.out.picard_alignmetrics, by: [0])
                      .join(Common_RNAseq_WF.out.picard_MD, by: [0])
                      .join(Common_RNAseq_WF.out.flagstat, by: [0])
                      .join(Common_RNAseq_WF.out.fastq_screen, by: [0])


multiqc_exome_input = Exome_common_WF.out.Fastqc_out
            .join(Exome_common_WF.out.verifybamid)
            .join(Exome_common_WF.out.flagstat)
            .join(Exome_common_WF.out.exome_final_bam)
            .join(Exome_common_WF.out.hsmetrics)
            .join(Exome_common_WF.out.krona)
            .join(Exome_common_WF.out.kraken)
            .join(Exome_common_WF.out.exome_qc)
            .join(Exome_common_WF.out.markdup_txt)

Combined_multiqc_input = multiqc_exome_input.merge(multiqc_rnaseq_input) { item1, item2 ->
    if (item1[0].id == item2[0].id && item1[0].casename == item2[0].casename) {
        return [[id: item1[0].id, casename: item1[0].casename]] + [item1[1..-1] + item2[1..-1]]
    } else {
        return null
    }

}
Multiqc(Combined_multiqc_input)

ch_versions = ch_versions.mix(Multiqc.out.versions)

combine_versions  = ch_versions.unique().collectFile(name: 'collated_versions.yml')

custom_versions_input = Multiqc.out.multiqc_report
        .combine(combine_versions).map{ meta, multiqc, version -> [meta, version] }
        .combine(Pipeline_version)

CUSTOM_DUMPSOFTWAREVERSIONS(custom_versions_input)


}
