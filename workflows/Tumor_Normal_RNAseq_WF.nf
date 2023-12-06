include {Common_RNAseq_WF} from './Common_RNAseq_WF'
include {Exome_common_WF} from './Exome_common_WF.nf'
include {MakeHotSpotDB} from '../modules/qc/plots'
include {Manta_Strelka} from '../subworkflows/Manta_Strelka.nf'
include {Mutect_WF} from '../subworkflows/Mutect.nf'
include {Exome_QC
        QC_summary_Patientlevel
        Multiqc} from '../modules/qc/qc.nf'
include {SnpEff
        Vcf2txt} from '../modules/misc/snpEff'
include {FormatInput} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {AddAnnotation_TN
        AddAnnotation_somatic_variants
        AddAnnotationFull_somatic_variants} from '../modules/annotation/annot'
include {UnionSomaticCalls} from '../modules/misc/UnionSomaticCalls.nf'
include {MutationalSignature
        Cosmic3Signature} from '../modules/misc/MutationalSignature.nf'
include {MutationBurden} from '../modules/misc/MutationBurden.nf'
include {Sequenza_annotation} from '../subworkflows/Sequenza_annotation'
include {Annotation_somatic} from '../subworkflows/Actionable_somatic.nf'
include {Annotation_germline} from '../subworkflows/Actionable_germline.nf'
include {Combine_variants
        VEP} from '../modules/annotation/VEP.nf'
include {DBinput_multiples } from '../modules/misc/DBinput'
//include {DBinput_multiples as DBinput_germline} from '../modules/misc/DBinput'
include {CNVkitPaired} from '../modules/cnvkit/CNVkitPaired'
include {CNVkit_png} from '../modules/cnvkit/CNVkitPooled'
include {TcellExtrect} from '../modules/misc/TcellExtrect'
include {Split_vcf
        Pvacseq
        Merge_Pvacseq_vcf} from '../modules/neoantigens/Pvacseq.nf'
include {CUSTOM_DUMPSOFTWAREVERSIONS} from '../modules/nf-core/dumpsoftwareversions/main.nf'

workflow Tumor_Normal_RNAseq_WF {

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
samples = Channel.fromPath("Tumor_RNAseq_Normal.csv")
.splitCsv(header:true)
.filter { row -> row.type == "tumor_DNA" || row.type == "normal_DNA" || row.type == "tumor_RNA" }
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
        exome: it[0].type == "normal_DNA" || it[0].type == "tumor_DNA"
        rnaseq:  it[0].type == "tumor_RNA"
}

samples_branch.rnaseq|Common_RNAseq_WF
samples_branch.exome|Exome_common_WF

exome_pileup_Status = Exome_common_WF.out.pileup.branch{
    normal: it[0].type == "normal_DNA"
    tumor:  it[0].type == "tumor_DNA"
}

//All Germline samples pileup in  [meta.id, meta, file] format
pileup_samples_normal_to_cross = exome_pileup_Status.normal.map{ meta, pileup -> [ meta.id, meta, pileup ] }

//All Tumor samples pileup  in [meta.id, meta, file] format
pileup_samples_tumor_to_cross = exome_pileup_Status.tumor.map{ meta, pileup -> [ meta.id, meta, pileup ] }

//All Tumor samples pileup  in [meta.id, meta, file] format
pileup_samples_rnaseq_to_cross = Common_RNAseq_WF.out.pileup.map{ meta, pileup -> [ meta.id, meta, pileup ] }


//Use cross to combine normal with tumor samples
pileup_pair = pileup_samples_normal_to_cross.cross(pileup_samples_tumor_to_cross)
            .map { normal, tumor ->
                def meta = [:]

                meta.id         = tumor[1].id
                meta.normal_id  = normal[1].lib
                meta.normal_type = normal[1].type
                meta.casename        = normal[1].casename
                meta.lib   = tumor[1].lib
                meta.type = tumor[1].type
                meta.T_sc  = tumor[1].sc
                meta.N_sc  = normal[1].sc

                [meta.id, meta, normal[2], tumor[2] ]
            }
            .cross(pileup_samples_rnaseq_to_cross)
            .map { exome, rnaseq ->
                def meta = exome[1]
                meta.rna_type = rnaseq[1].type

                [meta, [exome[2], exome[3],rnaseq[2]] ]
            }

MakeHotSpotDB(pileup_pair)

//tag the bam channel for Tumor and normal
bam_target_ch = Exome_common_WF.out.exome_final_bam.join(Exome_common_WF.out.target_capture_ch,by:[0])


bam_variant_calling_status = bam_target_ch.branch{
    normal: it[0].type == "normal_DNA"
    tumor:  it[0].type == "tumor_DNA"
}

// All Germline samples
bam_variant_calling_normal_to_cross = bam_variant_calling_status.normal.map{ meta, bam, bai, bed -> [ meta.id, meta, bam, bai, bed ] }

 // All tumor samples
bam_variant_calling_pair_to_cross = bam_variant_calling_status.tumor.map{ meta, bam, bai, bed -> [ meta.id, meta, bam, bai, bed ] }


bam_variant_calling_pair = bam_variant_calling_normal_to_cross.cross(bam_variant_calling_pair_to_cross)
            .map { normal, tumor ->
                def meta = [:]

                meta.id         = tumor[1].id
                meta.normal_id  = normal[1].lib
                meta.normal_type = normal[1].type
                meta.casename        = normal[1].casename
                meta.lib   = tumor[1].lib
                meta.type = tumor[1].type
                meta.T_sc  = tumor[1].sc
                meta.N_sc  = normal[1].sc

                [ meta, normal[2], normal[3], tumor[2], tumor[3],tumor[4] ]
            }

CNVkitPaired(
    bam_variant_calling_pair
    .combine(cnv_ref_access)
    .combine(genome)
    .combine(genome_fai)
    .combine(genome_dict)
)

ch_versions = Exome_common_WF.out.ch_versions.mix(CNVkitPaired.out.versions)

CNVkit_png(CNVkitPaired.out.cnvkit_pdf)

Manta_Strelka(bam_variant_calling_pair)

ch_versions = ch_versions.mix(Manta_Strelka.out.ch_versions)


SnpEff(Manta_Strelka.out.strelka_indel_raw_vcf
               .combine(dbNSFP2_4)
               .combine(dbNSFP2_4_tbi)
               .combine(Biowulf_snpEff_config)
               .combine(strelka_indelch)
    )

Vcf2txt(SnpEff.out.raw_snpeff.combine(strelka_indelch))

Mutect_WF(bam_variant_calling_pair)

ch_versions = ch_versions.mix(Mutect_WF.out.versions)

somatic_snpeff_input_ch = Mutect_WF.out.mutect_snpeff_snv_vcf2txt
        .join(Vcf2txt.out,by:[0])
        .join(Manta_Strelka.out.strelka_snpeff_snv_vcf2txt,by:[0])


HC_snpeff_snv_vcftxt_status = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.branch{
    normal: it[0].type == "normal_DNA"
    tumor:  it[0].type == "tumor_DNA"
}

//All Germline samples vcftxt in  [meta.id, meta, file] format
HC_snpeff_snv_vcftxt_samples_normal_to_cross = HC_snpeff_snv_vcftxt_status.normal.map{ meta, snpeff_snv_vcftxt  -> [ meta.id, meta, snpeff_snv_vcftxt ] }

//All Tumor samples vcftxt  in [meta.id, meta, file] format
HC_snpeff_snv_vcftxt_samples_tumor_to_cross = HC_snpeff_snv_vcftxt_status.tumor.map{ meta, snpeff_snv_vcftxt -> [ meta.id, meta, snpeff_snv_vcftxt ] }

HC_snpeff_snv_vcftxt_samples_rna_to_cross = Common_RNAseq_WF.out.snpeff_vcf.map{ meta, snpeff_snv_vcftxt -> [ meta.id, meta, snpeff_snv_vcftxt ] }


//Use cross to combine normal with tumor samples
snpeff_snv_vcftxt = HC_snpeff_snv_vcftxt_samples_normal_to_cross.cross(HC_snpeff_snv_vcftxt_samples_tumor_to_cross)
            .map { normal, tumor ->
                def meta = [:]

                meta.id         = tumor[1].id
                meta.normal_id  = normal[1].lib
                meta.normal_type = normal[1].type
                meta.casename        = normal[1].casename
                meta.lib   = tumor[1].lib
                meta.type = tumor[1].type
                meta.T_sc  = tumor[1].sc
                meta.N_sc  = normal[1].sc

                [ meta, normal[2], tumor[2] ]
            }
            .join(somatic_snpeff_input_ch,by:[0])
            .map{meta, hc_normal, hc_tumor, mutect, indels, snvs -> [ meta.id, meta, hc_normal, hc_tumor, mutect, indels, snvs ] }
            .cross(HC_snpeff_snv_vcftxt_samples_rna_to_cross)
            .map { exome, rnaseq ->
                def meta = exome[1]
                meta.rna_type = rnaseq[1].type

                [meta, [exome[2], exome[3], rnaseq[2], exome[4],exome[5],exome[6]] ]
            }

HC_snpeff_snv_vcftxt = snpeff_snv_vcftxt.map{ meta, files -> [meta, [files[0], files[1], files[2]] ] }
HC_snpeff_snv_vcftxt.view()

format_input_ch = snpeff_snv_vcftxt.join(MakeHotSpotDB.out,by:[0])

FormatInput(format_input_ch)

Annotation(FormatInput.out)

}
