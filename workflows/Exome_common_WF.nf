
include {BWA_picard} from '../subworkflows/Bwa_picard_bam_processing'
include {Exome_GATK} from '../subworkflows/Exome_GATK'
include {QC_exome_bam} from '../subworkflows/QC_exome_bam'
include {HLA_calls_exome} from '../subworkflows/HLA_calls_exome'
include {Kraken} from '../modules/qc/qc'
include {Krona} from '../modules/qc/qc'

workflow Exome_common_WF {

kraken_bacteria = Channel.of(file(params.kraken_bacteria, checkIfExists:true))

take:
     samples_exome_ch

main:

//Initiate empty channel for versions
ch_versions = Channel.empty()


BWA_picard(samples_exome_ch)

Kraken(samples_exome_ch
    .combine(kraken_bacteria)
)
Krona(Kraken.out.kraken_out)

HLA_calls_exome(samples_exome_ch)

capture_ch = BWA_picard.out.picard_MD
    .map { tuple ->
        def meta = tuple[0]
        def bam = tuple[1]
        def bai = tuple[2]
        def target_file = ''

        if (meta.sc == 'clin.ex.v1') {
            target_file = params.clin_ex_v1
        } else if (meta.sc == 'seqcapez.hu.ex.v3') {
            target_file = params.seqcapez.hu.ex.v3
        } else if (meta.sc == 'polya') {
            target_file = params.polya_target
        } else if (meta.sc == 'ribozero') {
            target_file = params.ribozero_target
        } else if (meta.sc == 'SmartRNA') {
            target_file = params.smartrna_target
        }
        } else if (meta.sc == 'agilent.v7') {
            target_file = params.agilent.v7_target
        }

        return [meta,target_file]
    }



Exome_GATK(BWA_picard.out.picard_MD,
           capture_ch
)

//update design channel
design_ch =  Exome_GATK.out.GATK_Exome_bam
   .map { tuple ->
        def meta = tuple[0]
        def bam = tuple[1]
        def bai = tuple[2]
        def design_file = ''

        if (meta.sc == 'clin.ex.v1') {
            design_file = params.clin_ex_v1_design
        } else if (meta.sc == 'seqcapez.hu.ex.v3') {
            design_file = params.seqcapez.hu.ex.v3
        } else if (meta.sc == 'polya') {
            design_file = params.polya_target
        } else if (meta.sc == 'ribozero') {
            design_file = params.ribozero_target
        } else if (meta.sc == 'SmartRNA') {
            design_file = params.smartrna_target
        }

        return [meta,design_file]
    }


QC_exome_bam(
    Exome_GATK.out.GATK_Exome_bam,
    BWA_picard.out.bwa_bam,
    capture_ch,
    design_ch
)


emit:
coverageplot = QC_exome_bam.out.coverageplot
pileup = QC_exome_bam.out.hotspot_pileup
snpeff_vcf =Exome_GATK.out.SnpEff_vcf
exome_final_bam = Exome_GATK.out.GATK_Exome_bam
loh = QC_exome_bam.out.loh
target_capture_ch = capture_ch
HC_snpeff_snv_vcf2txt = Exome_GATK.out.HC_snpeff_snv_vcf2txt
hlaminer_exome  = HLA_calls_exome.out.hlaminer_exome
seq2hla_exome = HLA_calls_exome.out.seq2hla_exome
mergehla_exome = HLA_calls_exome.out.mergehla_exome
exome_qc = QC_exome_bam.out.Exome_qc
}
