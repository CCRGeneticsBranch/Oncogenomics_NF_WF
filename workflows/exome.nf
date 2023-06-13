
include {BWA_picard} from '../subworkflows/Bwa_picard_bam_processing'
include {Exome_GATK} from '../subworkflows/Exome_GATK'
include {QC_exome_bam} from '../subworkflows/QC_exome_bam'

workflow EXOME {
  
//params.exome = "/data/khanlab/projects/Nextflow_dev/testing/exome_samplesheet.csv"
samples_exome = Channel.fromPath(params.samplesheet1)
.splitCsv(header:true)
.filter { row -> row.type == "Tumor" || row.type == "Normal" }
.map { row ->
    def meta = [:]
    meta.id    =  row.sample
    meta.lib   =  row.library
    meta.sc    =  row.sample_captures
    meta.casename  = row.casename 
    meta.type     = row.type
    def fastq_meta = []
    fastq_meta = [ meta,  file(row.read1), file(row.read2)  ]

    return fastq_meta
}
//samples_exome.view()
BWA_picard(samples_exome)

Exome_GATK(BWA_picard.out.picard_MD)

capture_ch = Exome_GATK.out.GATK_Exome_bam
    .map { tuple ->
        def meta = tuple[0]
        def bam = tuple[1]
        def bai = tuple[2]
        def target_file = ''
        
        if (meta.sc == 'clin.ex.v1') {
            target_file = params.clin_ex_v1
        } else if (meta.sc == 'polya_stranded') {
            target_file = params.polya_stranded_target
        } else if (meta.sc == 'polya') {
            target_file = params.polya_target
        } else if (meta.sc == 'ribozero') {
            target_file = params.ribozero_target
        } else if (meta.sc == 'SmartRNA') {
            target_file = params.smartrna_target
        }

        return [meta,target_file]
    }

//capture_ch.view()

QC_exome_bam(

    Exome_GATK.out.GATK_Exome_bam,
    BWA_picard.out.bwa_bam,
    capture_ch
)

}