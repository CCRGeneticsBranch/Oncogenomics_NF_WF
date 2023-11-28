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
include {Combine_customRNAQC
        RNAqc_TrancriptCoverage} from '../modules/qc/picard'
include {CircosPlot
        Genotyping_Sample} from '../modules/qc/qc'
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
Pipeline_version = Channel.from(params.Pipeline_version)


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
merge_fusion_anno_input = Fusion_Annotation.out.map{ meta, fusion -> [meta, [fusion]] }

Merge_fusion_annotation(merge_fusion_anno_input.combine(genome_version))

Combine_customRNAQC(Common_RNAseq_WF.out.rnalib_custum_qc.map{ meta, qc -> [meta, [qc]] })

RNAqc_TrancriptCoverage(Common_RNAseq_WF.out.picard_rnaseqmetrics.map{ meta, qc -> [meta, [qc]] })

MakeHotSpotDB(Common_RNAseq_WF.out.pileup.map{ meta, pileup -> [meta, [pileup]] })

//Run circos plot at case level
CircosPlot(Common_RNAseq_WF.out.loh.map{ meta, loh -> [meta, [loh]] })

Genotyping_Sample(Common_RNAseq_WF.out.gt.map{ meta, gt -> [meta, [gt]] })

formatinput_input_ch = Common_RNAseq_WF.out.snpeff_vcf.map{ meta, vcf -> [meta, [vcf]] }.join(MakeHotSpotDB.out)
FormatInput(formatinput_input_ch)


//Run Annotation subworkflow
Annotation(FormatInput.out)

ch_versions = Common_RNAseq_WF.out.ch_versions.mix(Annotation.out.version)

merged_ch = Common_RNAseq_WF.out.snpeff_vcf.join(Annotation.out.rare_annotation,by:[0])
AddAnnotation(merged_ch)

dbinput_anno_ch = AddAnnotation.out.map{ meta, txt -> [meta, [txt]] }
dbinput_snpeff_ch = Common_RNAseq_WF.out.snpeff_vcf.map{ meta, txt -> [meta, [txt]] }
dbinput_ch = dbinput_anno_ch.join(dbinput_snpeff_ch,by:[0])
DBinput(dbinput_ch)



multiqc_input = Common_RNAseq_WF.out.Fastqc_out.join(Common_RNAseq_WF.out.pileup, by: [0])
                      .join(Common_RNAseq_WF.out.coverageplot, by: [0])
                      .join(Common_RNAseq_WF.out.chimeric_junction, by: [0])
                      .join(Common_RNAseq_WF.out.rsem_genes, by: [0]).join(Common_RNAseq_WF.out.rnaseqc, by: [0])
                      .join(Common_RNAseq_WF.out.circos_plot, by: [0])
                      .join(Common_RNAseq_WF.out.strandedness, by: [0])


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
