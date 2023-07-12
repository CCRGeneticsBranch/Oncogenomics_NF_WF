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
include {Combine_variants} from '../modules/annotation/VEP.nf'
include {VEP} from '../modules/annotation/VEP.nf'

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

somatic_variants = Mutect_WF.out.mutect_raw_vcf
   .combine(Manta_Strelka.out.strelka_indel_raw_vcf,by:[0])
   .combine(Manta_Strelka.out.strelka_snvs_raw_vcf,by:[0])

HLA_normals = Exome_common_WF.out.hlaminer_exome.branch {Normal: it[0].type == "Normal"}
               .combine(Exome_common_WF.out.seq2hla_exome.branch {Normal: it[0].type == "Normal"},by:[0])
Combine_variants(
    somatic_variants,
    HLA_normals
)
//VEP(Combine_variants.out.combined_vcf_tmp)

}
