include {Exome_common_WF} from './Exome_common_WF.nf'
include {MakeHotSpotDB_TN} from '../modules/qc/plots'
include {Manta_Strelka} from '../subworkflows/Manta_Strelka.nf'
include {Mutect_WF} from '../subworkflows/Mutect.nf'
include {Exome_QC} from '../modules/qc/qc.nf'
include {SnpEff} from '../modules/misc/snpEff'
include {Vcf2txt} from '../modules/misc/snpEff'
include {FormatInput_TN} from '../modules/annotation/annot'
include {Annotation} from '../subworkflows/Annotation'
include {AddAnnotation_TN
        AddAnnotation_somatic_variants
        AddAnnotationFull_somatic_variants} from '../modules/annotation/annot'
include {UnionSomaticCalls} from '../modules/misc/UnionSomaticCalls.nf'
include {MutationalSignature} from '../modules/misc/MutationalSignature.nf'
include {Cosmic3Signature} from '../modules/misc/MutationalSignature.nf'
include {MutationBurden} from '../modules/misc/MutationBurden.nf'
include {Sequenza_annotation} from '../subworkflows/Sequenza_annotation'
include {Annotation_somatic} from '../subworkflows/Actionable_somatic.nf'
include {Annotation_germline} from '../subworkflows/Actionable_germline.nf'
include {Combine_variants} from '../modules/annotation/VEP.nf'
include {VEP} from '../modules/annotation/VEP.nf'
include {DBinput_multiples } from '../modules/misc/DBinput'
//include {DBinput_multiples as DBinput_germline} from '../modules/misc/DBinput'
include {QC_summary_Patientlevel} from '../modules/qc/qc'
include {CNVkitPaired} from '../modules/cnvkit/CNVkitPaired'
include {CNVkit_png} from '../modules/cnvkit/CNVkitPooled'
include {TcellExtrect} from '../modules/misc/TcellExtrect'
include {Split_vcf} from '../modules/neoantigens/Pvacseq.nf'
include {Pvacseq} from '../modules/neoantigens/Pvacseq.nf'
include {Merge_Pvacseq_vcf} from '../modules/neoantigens/Pvacseq.nf'
include {Multiqc} from '../modules/qc/qc'
include {CUSTOM_DUMPSOFTWAREVERSIONS} from '../modules/nf-core/dumpsoftwareversions/main.nf'


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
.filter { row -> row.type == "tumor_DNA" || row.type == "normal_DNA" }
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


pileup_Status = Exome_common_WF.out.pileup.branch{
    normal: it[0].type == "normal_DNA"
    tumor:  it[0].type == "tumor_DNA"
}

//All Germline samples pileup in  [meta.id, meta, file] format
pileup_samples_normal_to_cross = pileup_Status.normal.map{ meta, pileup -> [ meta.id, meta, pileup ] }

//All Tumor samples pileup  in [meta.id, meta, file] format
pileup_samples_tumor_to_cross = pileup_Status.tumor.map{ meta, pileup -> [ meta.id, meta, pileup ] }


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

                [ meta, normal[2], tumor[2] ]
            }

MakeHotSpotDB_TN(pileup_pair)


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


/*
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

//combined_HC_vcf_ch.view()
*/
HC_snpeff_snv_vcftxt_status = Exome_common_WF.out.HC_snpeff_snv_vcf2txt.branch{
    normal: it[0].type == "normal_DNA"
    tumor:  it[0].type == "tumor_DNA"
}

//All Germline samples pileup in  [meta.id, meta, file] format
HC_snpeff_snv_vcftxt_samples_normal_to_cross = HC_snpeff_snv_vcftxt_status.normal.map{ meta, snpeff_snv_vcftxt  -> [ meta.id, meta, snpeff_snv_vcftxt ] }

//All Tumor samples pileup  in [meta.id, meta, file] format
HC_snpeff_snv_vcftxt_samples_tumor_to_cross = HC_snpeff_snv_vcftxt_status.tumor.map{ meta, snpeff_snv_vcftxt -> [ meta.id, meta, snpeff_snv_vcftxt ] }


//Use cross to combine normal with tumor samples
HC_snpeff_snv_vcftxt_pair = HC_snpeff_snv_vcftxt_samples_normal_to_cross.cross(HC_snpeff_snv_vcftxt_samples_tumor_to_cross)
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


// Combine outputs from Mutect and strelka
somatic_snpeff_input_ch = Mutect_WF.out.mutect_snpeff_snv_vcf2txt
        .join(Vcf2txt.out,by:[0])
        .join(Manta_Strelka.out.strelka_snpeff_snv_vcf2txt,by:[0])


format_input_ch = somatic_snpeff_input_ch
        .join(HC_snpeff_snv_vcftxt_pair,by:[0])
        .join(MakeHotSpotDB_TN.out,by:[0])


FormatInput_TN(format_input_ch)

Annotation(FormatInput_TN.out)


addannotation_input_ch = HC_snpeff_snv_vcftxt_pair.join(Annotation.out.rare_annotation,by:[0])

AddAnnotation_TN(addannotation_input_ch)

addnnotation_somatic_variants_input_ch = somatic_snpeff_input_ch.join(Annotation.out.rare_annotation,by:[0])

AddAnnotation_somatic_variants(addnnotation_somatic_variants_input_ch)

addannotationfull_somatic_variants_input_ch = somatic_snpeff_input_ch.join(Annotation.out.final_annotation,by:[0])

AddAnnotationFull_somatic_variants(addannotationfull_somatic_variants_input_ch)

UnionSomaticCalls(AddAnnotationFull_somatic_variants.out)
//Test mutational signature only with full sample
//MutationalSignature(UnionSomaticCalls.out)



somatic_variants = Mutect_WF.out.mutect_raw_vcf
   .join(Manta_Strelka.out.strelka_indel_raw_vcf,by:[0])
   .join(Manta_Strelka.out.strelka_snvs_raw_vcf,by:[0])

//Test cosmic signature only with full sample
/*
Cosmic3Signature(
    somatic_variants
    .combine(cosmic_indel_rda)
    .combine(cosmic_genome_rda)
    .combine(cosmic_dbs_rda)
)
*/
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
                def meta = [:]

                meta.id         = vcf[1].id
                meta.normal_id  = vcf[1].normal_id
                meta.normal_type = vcf[1].normal_type
                meta.casename        = vcf[1].casename
                meta.lib   = vcf[1].lib
                meta.type = vcf[1].type
                meta.T_sc  = vcf[1].T_sc
                meta.N_sc  = vcf[1].N_sc

                [ meta, hla[2] ]
            }


pvacseq_input = split_vcf_files.combine(hla_with_updated_meta_ch,by:[0])

Pvacseq(pvacseq_input)

combined_pvacseq = Pvacseq.out.pvacseq_output_ch.groupTuple().map { meta, files -> [ meta, [*files] ] }


Merge_Pvacseq_vcf(combined_pvacseq)

addannotation_TN_combined_ch = AddAnnotation_TN.out.Tumor_hc_anno_txt.join(AddAnnotation_TN.out.Normal_hc_anno_txt,by:[0])

somatic_group               = Channel.from("somatic")
germline_group               = Channel.from("germline")

dbinput_somatic = AddAnnotation_somatic_variants.out
            .join(somatic_snpeff_input_ch,by:[0])
            .join(HC_snpeff_snv_vcftxt_pair,by:[0])
            .join(addannotation_TN_combined_ch,by:[0])


DBinput_multiples(dbinput_somatic.combine(somatic_group).combine(germline_group))

tumor_target_capture = bam_variant_calling_pair.map {meta, nbam, nbai, tbam, tbai, bed -> [ meta, bed ] }


Sequenza_annotation(
    bam_variant_calling_pair,
    tumor_target_capture)


ch_versions = ch_versions.mix(Sequenza_annotation.out.versions)


tcellextrect_input = Exome_common_WF.out.exome_final_bam
                    .join(Exome_common_WF.out.target_capture_ch,by:[0])
                    .combine(Sequenza_annotation.out.alternate)
                    .combine(genome_version_tcellextrect)

TcellExtrect(tcellextrect_input)



highconfidence_somatic_threshold = tumor_target_capture
   .map {tuple ->
        def meta = tuple[0]
        def bed = tuple[1]
        def Normal = ''
        def Tumor = ''
        def VAF =  ''
        if (meta.T_sc == 'clin.ex.v1' || meta.T_sc == 'nextera.ex.v1'|| meta.T_sc == 'vcrome2.1_pkv2' || meta.T_sc == 'seqcapez.hu.ex.v3' || meta.T_sc == 'seqcapez.hu.ex.utr.v1' || meta.T_sc == 'agilent.v7'|| meta.T_sc == 'panel_paed_v5_w5.1') {
            Normal = params.highconfidence_somatic_threshold['threshold_1']['Normal']
            Tumor = params.highconfidence_somatic_threshold['threshold_1']['Tumor']
            VAF = params.highconfidence_somatic_threshold['threshold_1']['VAF']
        } else if (meta.T_sc == 'clin.snv.v1'|| meta.T_sc == 'clin.snv.v2') {
            Normal = params.highconfidence_somatic_threshold['threshold_2']['Normal']
            Tumor = params.highconfidence_somatic_threshold['threshold_2']['Tumor']
            VAF = params.highconfidence_somatic_threshold['threshold_2']['VAF']
        } else if (meta.T_sc == 'seqcapez.rms.v1') {
            Normal = params.highconfidence_somatic_threshold['threshold_3']['Normal']
            Tumor = params.highconfidence_somatic_threshold['threshold_3']['Tumor']
            VAF = params.highconfidence_somatic_threshold['threshold_3']['VAF']
        } else if (meta.T_sc == 'wholegenome'){
            Normal = params.highconfidence_somatic_threshold['threshold_4']['Normal']
            Tumor = params.highconfidence_somatic_threshold['threshold_4']['Tumor']
            VAF = params.highconfidence_somatic_threshold['threshold_4']['VAF']
        }
        return [meta,Normal,Tumor,VAF]
   }

targetbp_MB_ch = tumor_target_capture
    .map { tuple ->
        def meta = tuple[0]
        def bed = tuple[1]
        def targetbp_mb = ''

        if (meta.T_sc == 'clin.ex.v1') {
            targetbp_mb = params.clin_ex_v1_MB
        } else if (meta.T_sc == 'seqcapez.hu.ex.v3') {
            targetbp_mb = params.seqcapez.hu.ex.v3_MB
        } else if (meta.T_sc == 'agilent.v7') {
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

qc_summary_ch = combinelibraries(Exome_common_WF.out.exome_qc)

QC_summary_Patientlevel(qc_summary_ch)




multiqc_input = Exome_common_WF.out.Fastqc_out
            .join(Exome_common_WF.out.verifybamid)
            .join(Exome_common_WF.out.flagstat)
            .join(Exome_common_WF.out.exome_final_bam)
            .join(Exome_common_WF.out.hsmetrics)
            .join(Exome_common_WF.out.krona)
            .join(Exome_common_WF.out.kraken)
            .join(Exome_common_WF.out.exome_qc)
            .join(Exome_common_WF.out.markdup_txt)


multiqc_status = multiqc_input.branch{
    normal: it[0].type == "normal_DNA"
    tumor:  it[0].type == "tumor_DNA"
}


multiqc_channel = multiqc_status.tumor.merge(multiqc_status.normal) { item1, item2 ->
    // Custom logic to merge items based on common id
    if (item1[0].id == item2[0].id && item1[0].casename == item2[0].casename) {
        return [[id: item1[0].id, casename: item1[0].casename]] + [item1[1..-1] + item2[1..-1]]
    } else {
        return null // or handle differently if ids don't match
    }
}
//multiqc_channel.view()
Multiqc(multiqc_channel)



ch_versions = ch_versions.mix(Multiqc_TN.out.versions)
/*
CUSTOM_DUMPSOFTWAREVERSIONS (
        tumor_meta,
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
        )

*/
}
