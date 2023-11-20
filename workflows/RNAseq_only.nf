//import workflows
include {Common_RNAseq_WF} from './Common_RNAseq_WF'

//import modules
include {MakeHotSpotDB} from  '../modules/qc/plots'
include {FormatInput} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {Multiqc} from '../modules/qc/qc'
include {Allstepscomplete} from '../modules/misc/Allstepscomplete'
include {AddAnnotation} from '../modules/annotation/annot'
include {DBinput} from '../modules/misc/DBinput'
include {RNAqc_TrancriptCoverage} from '../modules/qc/picard'
include {CircosPlot} from '../modules/qc/qc'
include {Actionable_variants} from '../modules/Actionable.nf'
include {Actionable_fusion} from '../modules/Actionable.nf'
include {Fusion_Annotation} from '../modules/annotation/Fusion_Annotation'
include {Merge_fusion_annotation} from '../modules/annotation/Fusion_Annotation'
include {CUSTOM_DUMPSOFTWAREVERSIONS} from '../modules/nf-core/dumpsoftwareversions/main.nf'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow RNAseq_only {


//config files
combined_gene_list = Channel.of(file(params.combined_gene_list, checkIfExists:true))
somatic_actionable_sites = Channel.of(file(params.somatic_actionable_sites, checkIfExists:true))
group               = Channel.from("rnaseq")
pfamdb  = Channel.of(file(params.pfamdb, checkIfExists:true))
genome  = Channel.of(file(params.genome, checkIfExists:true))
genome_version_fusion_annotation =  Channel.from(params.genome_version_fusion_annotation)
genome_version = Channel.from(params.genome_version)


//create a sample channel using meta hashmap
samples_rnaseq = Channel.fromPath("RNAseq.csv")
.splitCsv(header:true)
.filter { row -> row.type == "tumor_RNA" || row.type == "cell_line_RNA"}
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

//Run Common RNAseq WF, this runs all the steps from Cutadapt to GATK at library level
Common_RNAseq_WF(samples_rnaseq)
actionable_fusion_input = Common_RNAseq_WF.out.fusion_calls.map{ meta, fusion -> [meta, [fusion]] }
Actionable_fusion(actionable_fusion_input)


Fusion_Annotation_input = Common_RNAseq_WF.out.rsem_isoforms
                        .join(Common_RNAseq_WF.out.fusion_calls, by:[0])
                        .combine(pfamdb)
                        .combine(genome)
                        .combine(genome_version_fusion_annotation)
                        .combine(genome_version)


Fusion_Annotation(Fusion_Annotation_input)
/*


Merge_fusion_annotation(
    Fusion_Annotation.out.map { tuple -> tuple[1] },
    Fusion_Annotation.out.map { tuple -> tuple[0] },
    genome_version
)


//Converting meta channel to list channel
rnalib_qc_list_ch = Common_RNAseq_WF.out.rnalib_custum_qc.map { tuple -> tuple[1] }
rnaseqmetrics_list_ch = Common_RNAseq_WF.out.picard_rnaseqmetrics.map { tuple -> tuple[1] }
rnaseqmetrics_meta_ch = Common_RNAseq_WF.out.picard_rnaseqmetrics.map { tuple -> tuple[0] }

//Run custom RNA QC
RNAqc_TrancriptCoverage(
           rnalib_qc_list_ch,
           rnaseqmetrics_list_ch,
           rnaseqmetrics_meta_ch
)

//gather channels for makehotspotdb
pileup_input_ch = Common_RNAseq_WF.out.pileup.map { tuple -> tuple[1] }
pileup_meta_ch =Common_RNAseq_WF.out.pileup.map { tuple -> tuple[0] }

//Run Makehotspotdb
MakeHotSpotDB(pileup_input_ch,
             pileup_meta_ch
)

//Run FormatInput
formatinput_snpeff_ch = Common_RNAseq_WF.out.snpeff_vcf.map { tuple -> tuple.drop(1) }
FormatInput(
    formatinput_snpeff_ch,
    MakeHotSpotDB.out
)

//gather channels for circosplot
merged_loh_ch = Common_RNAseq_WF.out.loh.map { tuple -> tuple[1] }
meta_merged_loh = Common_RNAseq_WF.out.loh.map { tuple -> tuple[0] }

//Run circos plot at case level
CircosPlot(
    merged_loh_ch,
    meta_merged_loh
)

//Run Annotation subworkflow
Annotation(FormatInput.out)

merged_ch = Common_RNAseq_WF.out.snpeff_vcf.combine(Annotation.out.rare_annotation,by:[0])
AddAnnotation(merged_ch)

dbinput_snpeff_ch = Common_RNAseq_WF.out.snpeff_vcf.map{ tuple -> tuple.drop(1) }
dbinput_annot_ch = AddAnnotation.out.map{ tuple -> tuple.drop(1) }
dbinput_meta_ch = AddAnnotation.out.map { tuple -> tuple[0] }

//Run DBinput
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

//Common_RNAseq_WF.out.rsem_counts.view()

multiqc_input = Common_RNAseq_WF.out.Fastqc_out.join(Common_RNAseq_WF.out.pileup, by: [0])
                      .join(Common_RNAseq_WF.out.coverageplot, by: [0])
                      .join(Common_RNAseq_WF.out.chimeric_junction, by: [0])
                      .join(Common_RNAseq_WF.out.rsem_genes, by: [0]).join(Common_RNAseq_WF.out.rnaseqc, by: [0])
                      .join(Common_RNAseq_WF.out.circos_plot, by: [0])
Multiqc(multiqc_input.map { tuple -> tuple.drop(1) },
        multiqc_input.map { tuple -> tuple[0] })

final_inputs = Common_RNAseq_WF.out.fusion_calls.join(Annotation.out.final_annotation,by:[0]).join(Multiqc.out,by:[0])
Allstepscomplete(final_inputs)

Common_RNAseq_WF.out.ch_versions.unique().collectFile(name: 'collated_versions.yml').view()


CUSTOM_DUMPSOFTWAREVERSIONS (
        Common_RNAseq_WF.out.ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )

*/
}
