include {Exome_common_WF} from './Exome_common_WF.nf'
include {MakeHotSpotDB} from '../modules/qc/plots'
include {Manta_Strelka} from '../subworkflows/Manta_Strelka.nf'
include {Mutect_WF} from '../subworkflows/Mutect.nf'
include {Exome_QC} from '../modules/qc/qc.nf'
include {SnpEff} from '../modules/misc/snpEff'
include {Vcf2txt} from '../modules/misc/snpEff'
include {FormatInput} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {AddAnnotation} from '../modules/annotation/annot'
include {AddAnnotation_somatic_variants} from '../modules/annotation/annot'
include {AddAnnotationFull_somatic_variants} from '../modules/annotation/annot'
include {UnionSomaticCalls} from '../modules/misc/UnionSomaticCalls.nf'
include {MutationalSignature} from '../modules/misc/MutationalSignature.nf'
include {Cosmic3Signature} from '../modules/misc/MutationalSignature.nf'
include {MutationBurden} from '../modules/misc/MutationBurden.nf'
include {Sequenza_annotation} from '../subworkflows/Sequenza_annotation'
include {Annotation_somatic} from '../subworkflows/Actionable_somatic.nf'
include {Annotation_germline} from '../subworkflows/Actionable_germline.nf'
include {Combine_variants} from '../modules/annotation/VEP.nf'
include {VEP} from '../modules/annotation/VEP.nf'
include {QC_summary_Patientlevel} from '../modules/qc/qc'
include {CNVkitPaired} from '../modules/cnvkit/CNVkitPaired'
include {CNVkit_png} from '../modules/cnvkit/CNVkitPooled'
include {TcellExtrect} from '../modules/misc/TcellExtrect'
include {Split_vcf} from '../modules/neoantigens/Pvacseq.nf'
include {Pvacseq} from '../modules/neoantigens/Pvacseq.nf'



def combinelibraries(inputData) {
    def processedData = inputData.map { meta, file ->
        meta2 = [
            id: meta.id,
            casename: meta.casename
        ]
        [meta2, file]
    }.groupTuple()
     .map { meta, files -> [meta, *files] }
     .filter { tuple ->
        tuple.size() > 2
     }
    return processedData
}

workflow Tumor_Normal_WF {

    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    dbsnp_138_b37_vcf       = Channel.of(file(params.dbsnp, checkIfExists:true))
    cosmic_v67_hg19_vcf     = Channel.of(file(params.cosmic_v67_hg19_vcf, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
    strelka_config          = Channel.of(file(params.strelka_config, checkIfExists:true))
    strelka_indelch         = Channel.from("strelka.indels")
    strelka_snvsch          = Channel.from("strelka.snvs")
    mutect_ch               = Channel.from("MuTect")
    dbNSFP2_4             = Channel.of(file(params.dbNSFP2_4, checkIfExists:true))
    dbNSFP2_4_tbi         = Channel.of(file(params.dbNSFP2_4_tbi, checkIfExists:true))
    Biowulf_snpEff_config  = Channel.of(file(params.Biowulf_snpEff_config, checkIfExists:true))
    vep_cache              = Channel.of(file(params.vep_cache, checkIfExists:true))
    cosmic_indel_rda       = Channel.of(file(params.cosmic_indel_rda, checkIfExists:true))
    cosmic_genome_rda      = Channel.of(file(params.cosmic_genome_rda, checkIfExists:true))
    cosmic_dbs_rda         = Channel.of(file(params.cosmic_dbs_rda, checkIfExists:true))
    cnv_ref_access         = Channel.of(file(params.cnv_ref_access, checkIfExists:true))
    genome_version_tcellextrect         = Channel.of(params.genome_version_tcellextrect)

// Parse the samplesheet to generate fastq tuples
samples_exome = Channel.fromPath("Tumor_Normal.csv")
.splitCsv(header:true)
.filter { row -> row.type == "Tumor" || row.type == "Normal" }
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

//Run the common exome workflow, this workflow runs all the steps from BWA to GATK and QC steps
Exome_common_WF(samples_exome)



//Create a combined channel of libraries pileup  to generate hotspotdb at Patient-case level
Exome_common_WF.out.pileup.map { meta, file ->
    meta2 = [
        id: meta.id,
        casename: meta.casename
    ]
    [ meta2, file ]
  }.groupTuple()
   .map { meta, files -> [ meta, *files ] }
   .filter { tuple ->
    tuple.size() > 2
  }
   .set { combined_pileup_ch }
 

pileup_input_ch = combined_pileup_ch.map { tuple -> tuple.drop(1) }  
pileup_meta_ch = combined_pileup_ch.map { tuple -> tuple[0] }

MakeHotSpotDB(pileup_input_ch,
                   pileup_meta_ch
)

//tag the bam channel for Tumor 
bam_target_ch = Exome_common_WF.out.exome_final_bam.combine(Exome_common_WF.out.target_capture_ch,by:[0])
tumor_bam_channel = bam_target_ch.branch { 
    Tumor: it[0].type == "Tumor"
    Normal: it[0].type == "Normal"
}

CNVkitPaired(
    tumor_bam_channel.Tumor,
    tumor_bam_channel.Normal.map { tuple -> tuple.take(tuple.size() - 1) },
    cnv_ref_access,
    genome,
    genome_fai,
    genome_dict
)
CNVkit_png(CNVkitPaired.out.cnvkit_pdf)

Manta_Strelka(
    tumor_bam_channel.Tumor,
    tumor_bam_channel.Normal
)

SnpEff(Manta_Strelka.out.strelka_indel_raw_vcf
               .combine(dbNSFP2_4)
               .combine(dbNSFP2_4_tbi)
               .combine(Biowulf_snpEff_config)
               .combine(strelka_indelch)
    )

Vcf2txt(SnpEff.out.combine(strelka_indelch))

Mutect_WF(
    tumor_bam_channel.Tumor,
    tumor_bam_channel.Normal
)

Exome_common_WF.out.HC_snpeff_snv_vcf2txt.map { meta, file ->
    meta2 = [
        id: meta.id,
        casename: meta.casename
    ]
    [ meta2, file ]
  }.groupTuple()
   .map { meta, files -> [ meta, *files ] }
   .filter { tuple ->
    tuple.size() > 2
  }
   .set { combined_HC_vcf_ch }

Format_input_ch =  combined_HC_vcf_ch.map { tuple -> tuple.drop(1) }
        .combine(Mutect_WF.out.mutect_snpeff_snv_vcf2txt.map { tuple -> tuple.drop(1) })
        .combine(Vcf2txt.out.map { tuple -> tuple.drop(1) })
        .combine(Manta_Strelka.out.strelka_snpeff_snv_vcf2txt.map { tuple -> tuple.drop(1) })

FormatInput(
        Format_input_ch,
        MakeHotSpotDB.out
)

Annotation(FormatInput.out)
AddAnnotation_input_ch = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.combine(Annotation.out.rare_annotation).map { tuple ->[tuple[0], tuple[1], tuple[3]]}
AddAnnotation(AddAnnotation_input_ch)

somatic_variants_txt =Mutect_WF.out.mutect_snpeff_snv_vcf2txt
                            .combine(Vcf2txt.out,by:[0])
                            .combine(Manta_Strelka.out.strelka_snpeff_snv_vcf2txt,by:[0])
                
AddAnnotation_somatic_variants(
    somatic_variants_txt,
    Annotation.out.rare_annotation
)

AddAnnotationFull_somatic_variants(
    somatic_variants_txt,
    Annotation.out.final_annotation
)

UnionSomaticCalls(AddAnnotationFull_somatic_variants.out)
//Test mutational signature only with full sample
MutationalSignature(UnionSomaticCalls.out) 


somatic_variants = Mutect_WF.out.mutect_raw_vcf
   .combine(Manta_Strelka.out.strelka_indel_raw_vcf,by:[0])
   .combine(Manta_Strelka.out.strelka_snvs_raw_vcf,by:[0])

Cosmic3Signature(
    somatic_variants,
    cosmic_indel_rda,
    cosmic_genome_rda,
    cosmic_dbs_rda
)
Combine_variants(
    somatic_variants,
    Exome_common_WF.out.mergehla_exome.branch {Normal: it[0].type == "Normal"}
)

VEP(Combine_variants.out.combined_vcf_tmp.combine(vep_cache))

Split_vcf(VEP.out)



dbinput_somatic_annot = AddAnnotation_somatic_variants.out.map{ tuple -> tuple.drop(1) }
dbinput_somatic_snpeff = somatic_variants_txt.map{ tuple -> tuple.drop(1) }
dbinput_HC_snpeff = combined_HC_vcf_ch.map{ tuple -> tuple.drop(1) }
dbinput_meta_normal = (AddAnnotation.out.branch {Normal: it[0].type == "Normal"}.map { tuple -> tuple[0] })
dbinput_meta_tumor = (AddAnnotation.out.branch {Tumor: it[0].type == "Tumor"}.map { tuple -> tuple[0] })
/*
Annotation_somatic(
   dbinput_somatic_annot,
   dbinput_somatic_snpeff,
   dbinput_HC_snpeff,
   dbinput_meta_tumor,
   dbinput_meta_normal,
   Annotation.out.rare_annotation

)
*/
AddAnnotation.out.map { meta, file ->
    meta2 = [
        id: meta.id,
        casename: meta.casename
    ]
    [ meta2, file ]
  }.groupTuple()
   .map { meta, files -> [ meta, *files ] }
   .filter { tuple ->
    tuple.size() > 2
  }
   .set { dbinput_HC_annot_ch }

dbinput_HC_annot_ch = dbinput_HC_annot_ch.map{ tuple -> tuple.drop(1) }
/*
Annotation_germline(
   dbinput_HC_annot_ch,
   dbinput_somatic_snpeff,
   dbinput_HC_snpeff,
   dbinput_meta_tumor,
   dbinput_meta_normal,
   Annotation.out.rare_annotation,
   Annotation_somatic.out.dbinput_somatic

)
*/
tumor_target_capture = Exome_common_WF.out.target_capture_ch.branch { Tumor: it[0].type == "Tumor"}

Sequenza_annotation(
    tumor_bam_channel.Tumor,
    tumor_bam_channel.Normal,
    tumor_target_capture
)

tcellextrect_input = Exome_common_WF.out.exome_final_bam.combine(Exome_common_WF.out.target_capture_ch,by:[0]).combine(Sequenza_annotation.out.alternate).combine(genome_version_tcellextrect)

TcellExtrect(tcellextrect_input)

highconfidence_somatic_threshold = tumor_target_capture
   .map {tuple ->
        def meta = tuple[0]
        def bam = tuple[1]
        def Normal = ''
        def Tumor = ''
        def VAF =  ''
        if (meta.sc == 'clin.ex.v1' || meta.sc == 'nextera.ex.v1'|| meta.sc == 'vcrome2.1_pkv2' || meta.sc == 'seqcapez.hu.ex.v3' || meta.sc == 'seqcapez.hu.ex.utr.v1' || meta.sc == 'agilent.v7'|| meta.sc == 'panel_paed_v5_w5.1') {
            Normal = params.highconfidence_somatic_threshold['threshold_1']['Normal']
            Tumor = params.highconfidence_somatic_threshold['threshold_1']['Tumor']
            VAF = params.highconfidence_somatic_threshold['threshold_1']['VAF']
        } else if (meta.sc == 'clin.snv.v1'|| meta.sc == 'clin.snv.v2') {
            Normal = params.highconfidence_somatic_threshold['threshold_2']['Normal']
            Tumor = params.highconfidence_somatic_threshold['threshold_2']['Tumor']
            VAF = params.highconfidence_somatic_threshold['threshold_2']['VAF']
        } else if (meta.sc == 'seqcapez.rms.v1') {
            Normal = params.highconfidence_somatic_threshold['threshold_3']['Normal']
            Tumor = params.highconfidence_somatic_threshold['threshold_3']['Tumor']
            VAF = params.highconfidence_somatic_threshold['threshold_3']['VAF']
        } else if (meta.sc == 'wholegenome'){
            Normal = params.highconfidence_somatic_threshold['threshold_4']['Normal']
            Tumor = params.highconfidence_somatic_threshold['threshold_4']['Tumor']
            VAF = params.highconfidence_somatic_threshold['threshold_4']['VAF']
        }
        return [meta,Normal,Tumor,VAF] 
   }

MutationBurden(
    AddAnnotationFull_somatic_variants.out,
    params.clin_ex_v1_MB,
    highconfidence_somatic_threshold,
    mutect_ch,
    strelka_indelch,
    strelka_snvsch
)


def qc_summary_ch = combinelibraries(Exome_common_WF.out.exome_qc)

QC_summary_Patientlevel(qc_summary_ch)


//data/khanlab2/NF_benchmarking/work.vg_NCI0439/37/53b3dd612c5d2ec23959eafa24fde8/test/split
}
