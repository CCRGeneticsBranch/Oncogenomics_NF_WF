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
include {Actionable_RNAseq} from '../modules/Actionable.nf'
include {Actionable_fusion} from '../modules/Actionable.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow RNAseq_only {


//config files
combined_gene_list = Channel.of(file(params.combined_gene_list, checkIfExists:true))
somatic_actionable_sites = Channel.of(file(params.somatic_actionable_sites, checkIfExists:true))

//create a sample channel using meta hashmap
samples_rnaseq = Channel.fromPath("RNAseq.csv")
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

//Run Common RNAseq WF, this runs all the steps from Cutadapt to GATK at library level
Common_RNAseq_WF(samples_rnaseq)

//Create actionable fusions
Actionable_fusion(
Common_RNAseq_WF.out.fusion_calls.map { tuple -> tuple[1] },
Common_RNAseq_WF.out.fusion_calls.map { tuple -> tuple[0] }
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

Actionable_RNAseq(DBinput.out
       .combine(Annotation.out.rare_annotation,by:[0])
       .combine(combined_gene_list)
       .combine(somatic_actionable_sites)
)

//Common_RNAseq_WF.out.rsem_counts.view()

multiqc_input = Common_RNAseq_WF.out.Fastqc_out.join(Common_RNAseq_WF.out.pileup, by: [0])
                      .join(Common_RNAseq_WF.out.coverageplot, by: [0])
                      .join(Common_RNAseq_WF.out.chimeric_junction, by: [0])
                      .join(Common_RNAseq_WF.out.rsem_genes, by: [0]).join(Common_RNAseq_WF.out.rnaseqc, by: [0])
                      .join(Common_RNAseq_WF.out.circos_plot, by: [0])
Multiqc(multiqc_input)

final_inputs = Common_RNAseq_WF.out.fusion_calls.join(Annotation.out.final_annotation,by:[0]).join(Multiqc.out,by:[0])
Allstepscomplete(final_inputs)


}