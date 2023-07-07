include {Exome_common_WF} from './Exome_common_WF.nf'
include {MakeHotSpotDB} from '../modules/qc/plots'
include {Manta_Strelka} from '../subworkflows/Manta_Strelka.nf'
include {Mutect_WF} from '../subworkflows/Mutect.nf'
include {Exome_QC} from '../modules/qc/qc.nf'
include {SnpEff} from '../modules/misc/snpEff'
include {Vcf2txt} from '../modules/misc/snpEff'
include {FormatInput} from '../modules/annotation/annot'

workflow Tumor_Normal_WF {

    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    dbsnp_138_b37_vcf       = Channel.of(file(params.dbsnp, checkIfExists:true))
    cosmic_v67_hg19_vcf     = Channel.of(file(params.cosmic_v67_hg19_vcf, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
    strelka_config          = Channel.of(file(params.strelka_config, checkIfExists:true))
    strelka_indelch         = Channel.from("strelka_indels")
    dbNSFP2_4             = Channel.of(file(params.dbNSFP2_4, checkIfExists:true))
    dbNSFP2_4_tbi         = Channel.of(file(params.dbNSFP2_4_tbi, checkIfExists:true))
    Biowulf_snpEff_config  = Channel.of(file(params.Biowulf_snpEff_config, checkIfExists:true))


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
//samples_exome.view()

Exome_common_WF(samples_exome)

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


bam_target_ch = Exome_common_WF.out.exome_final_bam.combine(Exome_common_WF.out.target_capture_ch,by:[0])


tumor_bam_channel = bam_target_ch.branch { 
    Tumor: it[0].type == "Tumor"
    Normal: it[0].type == "Normal"
}

Manta_Strelka(
    tumor_bam_channel.Tumor,
    tumor_bam_channel.Normal
)

SnpEff(Manta_Strelka.out.strelka_indel_vcf
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

}
