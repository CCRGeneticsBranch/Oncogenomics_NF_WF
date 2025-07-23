//import workflows
include {Common_RNAseq_WF} from './Common_RNAseq_WF'

//import modules
include {RNAqc_TrancriptCoverage
        Combine_customRNAQC} from '../modules/qc/picard'
include {MakeHotSpotDB
        Hotspot_Boxplot
        CoveragePlot} from '../modules/qc/plots'
include {FormatInput} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {AddAnnotation} from '../modules/annotation/annot'
include {Fusion_Annotation
        Merge_fusion_annotation} from '../modules/annotation/Fusion_Annotation'
include {CircosPlot
        Genotyping_Sample} from '../modules/qc/qc'
include {Actionable_fusion} from '../modules/Actionable.nf'
include {Actionable_variants} from '../modules/Actionable.nf'
include {DBinput_multiple_new} from '../modules/misc/DBinput'
include {Multiqc} from '../modules/qc/qc'
include {CUSTOM_DUMPSOFTWAREVERSIONS} from '../modules/nf-core/dumpsoftwareversions/main.nf'
include {Allstepscomplete} from '../modules/misc/Allstepscomplete'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



def combinelibraries(inputData) {
    def processedData = inputData.map { meta, file ->
        meta2 = [
            id: meta.id,
            casename: meta.casename,
            diagnosis: meta.diagnosis
        ]
        [meta2, file]
    }.groupTuple()
     .map { meta, files -> [meta, [*files]] }
}



workflow RNAseq_multiple_libs {

//config files
combined_gene_list = Channel.of(file(params.combined_gene_list, checkIfExists:true))
somatic_actionable_sites = Channel.of(file(params.somatic_actionable_sites, checkIfExists:true))
genome_version_fusion_annotation =  Channel.from(params.genome_version_fusion_annotation)
genome_version = Channel.from(params.genome_version)
pfamdb  = Channel.of(file(params.pfamdb, checkIfExists:true))
genome  = Channel.of(file(params.genome, checkIfExists:true))
group               = Channel.from("rnaseq")
genome_version_fusion_annotation =  Channel.from(params.genome_version_fusion_annotation)
genome_version = Channel.from(params.genome_version)
Pipeline_version = Channel.from(params.Pipeline_version)

take:
    multiple_rna_samplesheet

main:

//create a sample channel using meta hashmap
samples_rnaseq = multiple_rna_samplesheet
.splitCsv(header:true)
.filter { row -> row.type == "tumor_RNA" || row.type == "cell_line_RNA" }
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
ch_allcomplete = Channel.empty()
//Run Common RNAseq WF, this runs all the steps from Cutadapt to GATK at library level
Common_RNAseq_WF(samples_rnaseq)


Fusion_Annotation_input = Common_RNAseq_WF.out.rsem_isoforms
                        .join(Common_RNAseq_WF.out.fusion_calls, by:[0])
                        .combine(pfamdb)
                        .combine(genome)
                        .combine(genome_version_fusion_annotation)
                        .combine(genome_version)


Fusion_Annotation(Fusion_Annotation_input)

merge_fusion_anno_input = Fusion_Annotation.out.map { meta, file ->
                    new_meta = [
                    id: meta.id,
                    casename: meta.casename
                    ]
                return [ new_meta, file ]
                }
                .groupTuple()
                .map { meta, files -> [ meta, [*files ]] }


Merge_fusion_annotation(merge_fusion_anno_input.combine(genome_version))
ch_allcomplete = ch_allcomplete.mix( Merge_fusion_annotation.out.map { all -> all[1..-1] }.flatten())

actionable_fusion_input = combinelibraries(Common_RNAseq_WF.out.fusion_calls)

Actionable_fusion(actionable_fusion_input)
ch_allcomplete = ch_allcomplete.mix( Actionable_fusion.out.map { meta, file -> file } )

combine_customRNAQC_input = combinelibraries(Common_RNAseq_WF.out.rnalib_custum_qc)
Combine_customRNAQC(combine_customRNAQC_input)
ch_allcomplete = ch_allcomplete.mix( Combine_customRNAQC.out.map { meta, file -> file } )


RNAqc_TrancriptCoverage_input = combinelibraries(Common_RNAseq_WF.out.picard_rnaseqmetrics)
RNAqc_TrancriptCoverage(RNAqc_TrancriptCoverage_input)
ch_allcomplete = ch_allcomplete.mix( RNAqc_TrancriptCoverage.out.map { meta, file -> file } )


genotyping_input = combinelibraries(Common_RNAseq_WF.out.gt)
Genotyping_Sample(genotyping_input,
                Pipeline_version)
ch_allcomplete = ch_allcomplete.mix( Genotyping_Sample.out.map { all -> all[1..-1] }.flatten())


circos_input = combinelibraries(Common_RNAseq_WF.out.loh)
CircosPlot(circos_input)
ch_allcomplete = ch_allcomplete.mix( CircosPlot.out.map { meta, file -> file } )


hotspot_depth_input = combinelibraries(Common_RNAseq_WF.out.hotspot_depth)
Hotspot_Boxplot(hotspot_depth_input)
ch_allcomplete = ch_allcomplete.mix( Hotspot_Boxplot.out.map { meta, file -> file } )


coverage_plot_input = combinelibraries(Common_RNAseq_WF.out.coverage)
CoveragePlot(coverage_plot_input)
ch_allcomplete = ch_allcomplete.mix( CoveragePlot.out.map { meta, file -> file } )


makehotspotdb_input = combinelibraries(Common_RNAseq_WF.out.pileup)
MakeHotSpotDB(makehotspotdb_input)

combined_snpefflibs = combinelibraries(Common_RNAseq_WF.out.snpeff_vcf)

format_input = combined_snpefflibs.join(MakeHotSpotDB.out,by:[0])
FormatInput(format_input)

Annotation(FormatInput.out)

ch_versions = Common_RNAseq_WF.out.ch_versions.mix(Annotation.out.version)


add_annotation_input = Common_RNAseq_WF.out.snpeff_vcf.combine(Annotation.out.rare_annotation.map{meta, file -> [file]})

AddAnnotation(add_annotation_input)

annot_ch_dbinput = AddAnnotation.out.map{ meta, file -> [ meta.id, meta.casename, meta, file ] }
            .map { patient, casename, meta, file ->
            meta2 = [
                lib: meta.lib,
                sc: meta.sc,
                type: meta.type
            ]
            [patient,casename, meta2, file]
    }
    .reduce([[:], [], [], []]) { result, item ->
        def (patient, casename, meta, filePath) = item

        // Dynamically set the id from the meta data
        result[0] = [id: patient, casename: casename]

        // Append library and type to the result
        result[1].add(meta.lib)
        result[1].add(meta.type)

        // Append library and source class to the result
        result[2].add(meta.lib)
        result[2].add(meta.sc)

        // Append file paths
        result[3].add(filePath)

        return result
    }


DBinput_multiple_new(annot_ch_dbinput,
                    combined_snpefflibs.map{meta, file -> file})
ch_allcomplete = ch_allcomplete.mix( DBinput_multiple_new.out.map { meta, file -> file } )



multiqc_input = Common_RNAseq_WF.out.Fastqc_out.join(Common_RNAseq_WF.out.pileup, by: [0])
                   .join(Common_RNAseq_WF.out.chimeric_junction, by: [0])
                   .join(Common_RNAseq_WF.out.rsem_genes, by: [0])
                   .join(Common_RNAseq_WF.out.rnaseqc, by: [0])
                   .join(Common_RNAseq_WF.out.circos_plot, by: [0])
                   .join(Common_RNAseq_WF.out.strandedness, by: [0])
                   .join(Common_RNAseq_WF.out.rnalib_custum_qc, by: [0])
                   .join(Common_RNAseq_WF.out.picard_rnaseqmetrics, by: [0])
                   .join(Common_RNAseq_WF.out.picard_rnaseqmetrics_pdf, by: [0])
                   .join(Common_RNAseq_WF.out.picard_alignmetrics, by: [0])
                   .join(Common_RNAseq_WF.out.markdup_txt, by: [0])
                   .join(Common_RNAseq_WF.out.flagstat, by: [0])
                   .join(Common_RNAseq_WF.out.fastq_screen, by: [0])

multiqc_input.map { meta, fastqc, pileup, chimeric, rsem, rnaseqc, circos, strand, rna_qc, picardqc, picardqc_pdf, picardqc_metric, markdup, flagstat, fastq_screen ->
    meta2 = [
        id: meta.id,
        casename: meta.casename,
    ]
    [ meta2, fastqc, pileup, chimeric, rsem, rnaseqc, circos, strand, rna_qc, picardqc, picardqc_pdf, picardqc_metric, markdup, flagstat, fastq_screen ]
  }.groupTuple()
   .map { meta, fastqcs, pileups, chimerics, rsems, rnaseqcs, circoss, strands, rna_qcs, picardqcs, picardqc_pdfs, picardqc_metrics, markdups, flagstats, fastq_screens -> [ meta, [*fastqcs, *pileups, *chimerics, *rsems, *rnaseqcs, *circoss, *strands, *rna_qcs, *picardqcs, *picardqc_pdfs, *picardqc_metrics, *markdups, *flagstats, *fastq_screens ]] }
   .set { multiqc_ch }

Multiqc(multiqc_ch)
ch_versions = ch_versions.mix(Multiqc.out.versions)

combine_versions  = ch_versions.unique().collectFile(name: 'collated_versions.yml')

custom_versions_input = Multiqc.out.multiqc_report
        .combine(combine_versions).map{ meta, multiqc, version -> [meta, version] }
        .combine(Pipeline_version)

CUSTOM_DUMPSOFTWAREVERSIONS(custom_versions_input)

Allstepscomplete(CUSTOM_DUMPSOFTWAREVERSIONS.out.config,
                ch_allcomplete)

}
