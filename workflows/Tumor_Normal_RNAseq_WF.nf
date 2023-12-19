include {Common_RNAseq_WF} from './Common_RNAseq_WF'
include {Exome_common_WF} from './Exome_common_WF.nf'
include {Manta_Strelka} from '../subworkflows/Manta_Strelka.nf'
include {Mutect_WF} from '../subworkflows/Mutect.nf'
include {Exome_QC
        QC_summary_Patientlevel
        Multiqc} from '../modules/qc/qc.nf'
include {SnpEff
        Vcf2txt} from '../modules/misc/snpEff'
include {MakeHotSpotDB} from '../modules/qc/plots'
include {FormatInput} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {AddAnnotation_TN
        AddAnnotation_somatic_variants
        AddAnnotationFull_somatic_variants
        Expressed} from '../modules/annotation/annot'
include {UnionSomaticCalls} from '../modules/misc/UnionSomaticCalls.nf'
include {MutationalSignature
        Cosmic3Signature} from '../modules/misc/MutationalSignature.nf'
include {MutationBurden} from '../modules/misc/MutationBurden.nf'
include {Sequenza_annotation} from '../subworkflows/Sequenza_annotation'
include {Annotation_somatic} from '../subworkflows/Actionable_somatic.nf'
include {Annotation_germline} from '../subworkflows/Actionable_germline.nf'
include {Combine_variants
        VEP} from '../modules/annotation/VEP.nf'
include {DBinput} from '../modules/misc/DBinput'
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
                def meta = tumor[1]
                [
                    meta + [
                        N_sc: normal[1].sc,
                        normal_id: normal[1].lib,
                        normal_type: normal[1].type
        ],
        normal[2], normal[3], tumor[2], tumor[3],tumor[4]
                ]

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
HC_snpeff_snv_vcftxt_exome = HC_snpeff_snv_vcftxt_samples_normal_to_cross.cross(HC_snpeff_snv_vcftxt_samples_tumor_to_cross)
            .map { normal, tumor ->
                def meta = tumor[1]
                [
                    meta + [
                        N_sc: normal[1].sc,
                        normal_id: normal[1].lib,
                        normal_type: normal[1].type
        ],

        normal[2], tumor[2] ]
            }

HC_snpeff_somatic_snpeff_input_ch_to_cross = HC_snpeff_snv_vcftxt_exome.join(somatic_snpeff_input_ch,by:[0])
            .map{meta, hc_normal, hc_tumor, mutect, indels, snvs -> [ meta.id, meta, hc_normal, hc_tumor, mutect, indels, snvs ] }

snpeff_snv_vcftxt = HC_snpeff_somatic_snpeff_input_ch_to_cross.cross(HC_snpeff_snv_vcftxt_samples_rna_to_cross)
            .map { exome, rnaseq ->
                def meta = exome[1]
                [
                    meta + [
                rna_lib: rnaseq[1].lib,
                rna_type: rnaseq[1].type,
                RNA_sc: rnaseq[1].sc
        ],

        [exome[2], exome[3], rnaseq[2], exome[4],exome[5],exome[6]] ]
            }

HC_snpeff_snv_vcftxt = snpeff_snv_vcftxt.map{ meta, files -> [meta, [files[0], files[1], files[2]] ] }


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
pileup_pair_exome = pileup_samples_normal_to_cross.cross(pileup_samples_tumor_to_cross)
            .map { normal, tumor ->
                def meta = tumor[1]
                [ meta.id,
                    meta + [
                        N_sc: normal[1].sc,
                        normal_id: normal[1].lib,
                        normal_type: normal[1].type
        ],


        normal[2], tumor[2] ]
            }
pileup_pair =  pileup_pair_exome.cross(pileup_samples_rnaseq_to_cross)
            .map { exome, rnaseq ->
                def meta = exome[1]
                [
                    meta + [
                rna_lib: rnaseq[1].lib,
                rna_type: rnaseq[1].type,
                RNA_sc: rnaseq[1].sc
        ],

        [exome[2], exome[3],rnaseq[2]] ]
            }



MakeHotSpotDB(pileup_pair)

format_input_ch = snpeff_snv_vcftxt.join(MakeHotSpotDB.out,by:[0])


FormatInput(format_input_ch)


Annotation(FormatInput.out)

addannotation_input_ch = HC_snpeff_snv_vcftxt.join(Annotation.out.rare_annotation,by:[0])

AddAnnotation_TN(addannotation_input_ch)


//addnnotation_somatic_variants_input_ch = somatic_snpeff_input_ch.join(Annotation.out.rare_annotation,by:[0])
somatic_snpeff_input_ch_to_cross = somatic_snpeff_input_ch.map{meta, mutect, indels, snv  -> [ meta.id, meta, mutect, indels, snv ] }
Annotation_rare_to_cross = Annotation.out.rare_annotation.map{meta, rare -> [meta.id, meta, rare ] }

addnnotation_somatic_variants_input_ch = somatic_snpeff_input_ch_to_cross.cross(Annotation_rare_to_cross)
                        .map {snpeff, annotation ->
                        def meta = annotation[1]
                        [meta, snpeff[2],snpeff[3],snpeff[4],annotation[2]]
                        }

AddAnnotation_somatic_variants(addnnotation_somatic_variants_input_ch)


RNA_expressed_ch_to_cross = Common_RNAseq_WF.out.Bam.join(Common_RNAseq_WF.out.snpeff_vcf,by:[0])
                                .map{meta, bam, txt -> [meta.id, meta, bam, txt ] }

AddAnnotation_somatic_variants_expressed_ch_to_cross = AddAnnotation_somatic_variants.out.map{meta, mutect, indels, snv -> [ meta.id, meta, mutect, indels, snv ] }


expressed_input_ch = AddAnnotation_somatic_variants_expressed_ch_to_cross.cross(RNA_expressed_ch_to_cross)
                        .map {annotation, rnaseq ->
                        def meta = annotation[1]
                        [meta, annotation[2], annotation[3], annotation[4],[rnaseq[2],rnaseq[3]]]
                        }

Expressed(expressed_input_ch)

Annotation_final_to_cross = Annotation.out.final_annotation.map{meta, rare -> [meta.id, meta, rare ] }

addannotationfull_somatic_variants_input_ch = somatic_snpeff_input_ch_to_cross.cross(Annotation_final_to_cross)
                        .map {snpeff, annotation ->
                        def meta = annotation[1]
                        [meta, snpeff[2],snpeff[3],snpeff[4],annotation[2]]
                        }


AddAnnotationFull_somatic_variants(addannotationfull_somatic_variants_input_ch)

/*
//Expressed_somatic_out = Expressed.out.map{ meta, mutect, indels, snv -> [meta,[mutect, indels, snv]]}

HC_annotated = AddAnnotation_TN.out.Normal_hc_anno_txt
                .join(AddAnnotation_TN.out.Tumor_hc_anno_txt,by:[0])
                .join(AddAnnotation_TN.out.RNA_hc_anno_txt,by:[0])

HC_annotated_Expressed = HC_annotated.join(Expressed.out,by:[0])
                        .map{meta, normal_anno, tumor_anno, rna_anno, mutect, indels, snv -> [meta,[normal_anno, tumor_anno, rna_anno, mutect, indels, snv] ] }

dbinput_ch = HC_annotated_Expressed.join(snpeff_snv_vcftxt,by:[0])

DBinput(dbinput_ch)
UnionSomaticCalls(AddAnnotationFull_somatic_variants.out)
//Test mutational signature only with full sample
//MutationalSignature(UnionSomaticCalls.out)



somatic_variants = Mutect_WF.out.mutect_raw_vcf
   .join(Manta_Strelka.out.strelka_indel_raw_vcf,by:[0])
   .join(Manta_Strelka.out.strelka_snvs_raw_vcf,by:[0])
/*
Cosmic3Signature(
    somatic_variants
    .combine(cosmic_indel_rda)
    .combine(cosmic_genome_rda)
    .combine(cosmic_dbs_rda)
)

Combine_variants(somatic_variants)

VEP(Combine_variants.out.combined_vcf_tmp.combine(vep_cache))

ch_versions = ch_versions.mix(VEP.out.versions)

Split_vcf(VEP.out.vep_out)


split_vcf_files = Split_vcf.out.flatMap { meta, files -> files.collect { [meta, it] } }


mergehla_status = Exome_common_WF.out.mergehla_exome.branch{
    normal: it[0].type == "normal_DNA"
    tumor:  it[0].type == "tumor_DNA"
}

//All Germline samples HLA in  [meta.id, meta, file] format
mergehla_samples_normal_to_cross = mergehla_status.normal.map{ meta, mergedcalls  -> [ meta.id, meta, mergedcalls ] }

vep_out_to_cross = VEP.out.vep_out.map{ meta, vcf -> [ meta.id, meta, vcf ] }

hla_with_updated_meta_ch = vep_out_to_cross.cross(mergehla_samples_normal_to_cross)
            .map { vcf, hla ->
                def meta = vcf[1]
                [ meta, hla[2] ]
            }

pvacseq_input = split_vcf_files.combine(hla_with_updated_meta_ch,by:[0])

Pvacseq(pvacseq_input)

combined_pvacseq = Pvacseq.out.pvacseq_output_ch.groupTuple().map { meta, files -> [ meta, [*files] ] }
Merge_Pvacseq_vcf(combined_pvacseq)

tumor_target_capture = bam_variant_calling_pair.map {meta, nbam, nbai, tbam, tbai, bed -> [ meta, bed ] }
/*
Sequenza_annotation(
    bam_variant_calling_pair,
    tumor_target_capture)


ch_versions = ch_versions.mix(Sequenza_annotation.out.versions)

tcellextrect_input = Exome_common_WF.out.exome_final_bam
                    .join(Exome_common_WF.out.target_capture_ch,by:[0])
                    .combine(genome_version_tcellextrect)

TcellExtrect(tcellextrect_input)
*/

/*
highconfidence_somatic_threshold = pileup_pair
   .map {tuple ->
        def meta = tuple[0]
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

targetbp_MB_ch = pileup_pair
    .map { tuple ->
        def meta = tuple[0]
        def targetbp_mb = ''

        if (meta.sc == 'clin.ex.v1') {
            targetbp_mb = params.clin_ex_v1_MB
        } else if (meta.sc == 'seqcapez.hu.ex.v3') {
            targetbp_mb = params.seqcapez.hu.ex.v3_MB
        } else if (meta.sc == 'agilent.v7') {
            targetbp_mb = params.agilent.v7_MB
        }

        return [meta,targetbp_mb]
    }

mutationburden_input_ch = AddAnnotationFull_somatic_variants.out
                    .join(targetbp_MB_ch,by:[0])
                    .join(highconfidence_somatic_threshold,by:[0])
                    .combine(mutect_ch)
                    .combine(strelka_indelch)
                    .combine(strelka_snvsch)

MutationBurden(mutationburden_input_ch)

*/

}
