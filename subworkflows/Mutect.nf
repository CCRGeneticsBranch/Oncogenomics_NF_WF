include {Mutect} from '../modules/Variant_analysis/Mutect'
include {Mutect2} from '../modules/Variant_analysis/Mutect'
include {Mutect_order} from '../modules/Variant_analysis/Mutect'
include {learn_read_orientation_model} from '../modules/Variant_analysis/Mutect'
include {filter_mutect} from '../modules/Variant_analysis/Mutect'
include {SnpEff} from '../modules/misc/snpEff'
include {Vcf2txt} from '../modules/misc/snpEff'

workflow Mutect_WF {

    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
    dbsnp_138_b37_vcf       = Channel.of(file(params.dbsnp, checkIfExists:true))
    dbNSFP2_4             = Channel.of(file(params.dbNSFP2_4, checkIfExists:true))
    dbNSFP2_4_tbi         = Channel.of(file(params.dbNSFP2_4_tbi, checkIfExists:true))
    Biowulf_snpEff_config  = Channel.of(file(params.Biowulf_snpEff_config, checkIfExists:true))
    mutect_ch               = Channel.from("MuTect")
take:

    bam_variant_calling_pair

main:

if (params.genome_v == "hg38") {

    // Filter to only paired samples that have a normal_id in meta
    bam_variant_calling_pair_with_normal = bam_variant_calling_pair
        .filter { it[0].normal_id }

    germline_resource     = Channel.of(file(params.gnomad, checkIfExists:true))
    germline_resource_idx = Channel.of(file(params.gnomad_index, checkIfExists:true))
    pon                   = Channel.of(file(params.tcga_pon, checkIfExists:true))
    pon_idx               = Channel.of(file(params.tcga_pon_index, checkIfExists:true))

    Mutect2(
        bam_variant_calling_pair_with_normal,
        genome,
        genome_fai,
        genome_dict,
        germline_resource,
        germline_resource_idx,
        pon,
        pon_idx
    )

    learn_read_orientation_model(Mutect2.out.f1r2)

    filter_mutect(
        Mutect2.out.mutect2_raw_vcf
            .join(Mutect2.out.mutect2_raw_vcf_tbi, by:[0])
            .join(Mutect2.out.mutect2_stats, by:[0])
            .join(learn_read_orientation_model.out.orientation_model, by:[0]),
        genome,
        genome_fai,
        genome_dict
    )

} else {

    cosmic_v67_hg19_vcf = Channel.of(file(params.cosmic_v67_hg19_vcf, checkIfExists:true))

    Mutect(
        bam_variant_calling_pair
        .combine(genome)
        .combine(genome_fai)
        .combine(genome_dict)
        .combine(dbsnp_138_b37_vcf)
        .combine(cosmic_v67_hg19_vcf)
    )

    Mutect_order(Mutect.out.mutect_raw_vcf)
}

    // Unified VCF channel from either filter_mutect (hg38) or Mutect_order (hg19)
    mutect_vcf_ch = params.genome_v == "hg38" ? filter_mutect.out.mutect2_filtered_vcf : Mutect_order.out

    SnpEff(mutect_vcf_ch
           .combine(dbNSFP2_4)
           .combine(dbNSFP2_4_tbi)
           .combine(Biowulf_snpEff_config)
           .combine(mutect_ch)
    )
    Vcf2txt(SnpEff.out.raw_snpeff.combine(mutect_ch))

emit:
    mutect_snpeff_snv_vcf2txt = Vcf2txt.out
    mutect_raw_vcf = mutect_vcf_ch
    mutect2_raw_vcf_tbi = params.genome_v == "hg38" ? filter_mutect.out.mutect2_filtered_vcf_tbi : Channel.empty()
    mutect2_stats = params.genome_v == "hg38" ? Mutect2.out.mutect2_stats : Channel.empty()
    f1r2 = params.genome_v == "hg38" ? Mutect2.out.f1r2 : Channel.empty()
    orientation_model = params.genome_v == "hg38" ? learn_read_orientation_model.out.orientation_model : Channel.empty()
    versions = params.genome_v == "hg38" ? Mutect2.out.versions : Mutect.out.versions

}
