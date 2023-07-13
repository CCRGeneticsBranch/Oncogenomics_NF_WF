//import workflows
include {Common_RNAseq_WF} from './Common_RNAseq_WF'

//import modules
include {RNAqc_TrancriptCoverage} from '../modules/qc/picard'
include {MakeHotSpotDB} from '../modules/qc/plots'
include {FormatInput} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {AddAnnotation} from '../modules/annotation/annot'
include {Fusion_Annotation} from '../modules/annotation/Fusion_Annotation'
include {CircosPlot} from '../modules/qc/qc'
include {Actionable_fusion} from '../modules/Actionable.nf'
include {Actionable_RNAseq} from '../modules/Actionable.nf'
include {DBinput_multiple} from '../modules/misc/DBinput'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow RNAseq_multiple_libs {

//config files
combined_gene_list = Channel.of(file(params.combined_gene_list, checkIfExists:true))
somatic_actionable_sites = Channel.of(file(params.somatic_actionable_sites, checkIfExists:true))


//create a sample channel using meta hashmap
samples_rnaseq = Channel.fromPath("RNA_lib.csv")
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

//create combined fusion calls channel of multiple libraries
Common_RNAseq_WF.out.fusion_calls.map { meta, file ->
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
   .set { combined_fusion_ch }

//Create actionable fusions
Actionable_fusion(
combined_fusion_ch.map { tuple -> tuple.drop(1) },
combined_fusion_ch.map { tuple -> tuple[0] }
)

//create combined isoform calls channel of multiple libraries
Common_RNAseq_WF.out.rsem_isoforms.map { meta, file ->
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
   .set { combined_rsem_ch }
/*
Fusion_Annotation(
    combined_rsem_ch.map { tuple -> tuple.drop(1) },
    Actionable_fusion.out
)
*/

//create combined qc channel of multiple libraries
Common_RNAseq_WF.out.rnalib_custum_qc.map { meta, file ->
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
   .set { combined_RNAqc_ch }

//create combined Picard qc channel of multiple libraries
Common_RNAseq_WF.out.picard_rnaseqmetrics.map { meta, file ->
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
   .set { combined_rnaseqmetrics_ch }

//gather channels for custum qc
combined_RNAqc_list_ch = combined_RNAqc_ch.map { tuple -> tuple.drop(1) }
combined_rnaseqmetrics_list_ch = combined_rnaseqmetrics_ch.map { tuple -> tuple.drop(1) }
combined_rnaseqmetrics_meta_ch = combined_rnaseqmetrics_ch.map { tuple -> tuple[0] }

//Custom RNA QC
RNAqc_TrancriptCoverage(
    combined_RNAqc_list_ch,
    combined_rnaseqmetrics_list_ch,
    combined_rnaseqmetrics_meta_ch
)

//create combined mpileup channel of multiple libraries
Common_RNAseq_WF.out.pileup.map { meta, file ->
    meta2 = [
        id: meta.id,
        casename: meta.casename,
        type: meta.type,
        sc: meta.sc
    ]
    [ meta2, file ]
  }.groupTuple()
   .map { meta, files -> [ meta, *files ] }
   .filter { tuple ->
    tuple.size() > 2
  }
   .set { combined_pileup_ch }

//gather channels for makehotspotdb
pileup_input_ch = combined_pileup_ch.map { tuple -> tuple.drop(1) }  
pileup_meta_ch = combined_pileup_ch.map { tuple -> tuple[0] }

//Run Makehotspotdb
MakeHotSpotDB(pileup_input_ch,
                   pileup_meta_ch
)

//create combined snpeff channel of multiple libraries
Common_RNAseq_WF.out.snpeff_vcf.map { meta, file ->
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
   .set { combined_snpeff_ch }


//Run FormatInput
FormatInput(
    combined_snpeff_ch.map { tuple -> tuple.drop(1) },
    MakeHotSpotDB.out
)

//Run Annotation subworkflow
Annotation(FormatInput.out)

//create combined loh channel of multiple libraries
Common_RNAseq_WF.out.loh.map { meta, file ->
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
   .set { combined_loh }

//Run circos plot at case level
CircosPlot(
    combined_loh.map { tuple -> tuple.drop(1) },
    combined_loh.map { tuple -> tuple[0] }
)

//create combined channel of snpeff.txt and rare_annotation
merged_ch = Common_RNAseq_WF.out.snpeff_vcf.combine(Annotation.out.rare_annotation)
updated_tuples = merged_ch.map { tuple ->
    [tuple[0], tuple[1], tuple[3]]
}
//Run AddAnnotation
AddAnnotation(updated_tuples)    


//create combined snpeff channel of multiple libraries
Common_RNAseq_WF.out.snpeff_vcf.map { meta, file ->
    meta2 = [
        id: meta.id,
        casename: meta.casename,
        type: meta.type,
        sc: meta.sc
    ]
    [ meta2, file ]
  }.groupTuple()
   .map { meta, files -> [ meta, *files ] }
   .filter { tuple ->
    tuple.size() > 2
  }
   .set { dbinput_combined_snpeff_txt }

//create combined Annotation channel of multiple libraries
AddAnnotation.out.map { meta, file ->
    meta2 = [
        id: meta.id,
        casename: meta.casename,
        type: meta.type,
        sc: meta.sc
    ]
    [ meta2, file ]
  }.groupTuple()
   .map { meta, files -> [ meta, *files ] }
   .filter { tuple ->
    tuple.size() > 2
  }
   .set { dbinput_combined_addannot_txt }

//Run DBinput
DBinput_multiple(
    dbinput_combined_addannot_txt,
    dbinput_combined_snpeff_txt
)

//Run Actionable_RNAseq to generate .rnaseq files
Actionable_RNAseq(DBinput_multiple.out
       .combine(Annotation.out.rare_annotation,by:[0])
       .combine(combined_gene_list)
       .combine(somatic_actionable_sites)
)

multiqc_input = Common_RNAseq_WF.out.Fastqc_out.join(Common_RNAseq_WF.out.pileup, by: [0])
                   .join(Common_RNAseq_WF.out.coverageplot, by: [0])
                   .join(Common_RNAseq_WF.out.chimeric_junction, by: [0])
                   .join(Common_RNAseq_WF.out.rsem_genes, by: [0]).join(Common_RNAseq_WF.out.rnaseqc, by: [0])
                   .join(Common_RNAseq_WF.out.circos_plot, by: [0])
//multiqc_input.view()


}

