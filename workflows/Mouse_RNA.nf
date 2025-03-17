/*
include {Cutadapt_se
         Star_se
         Rsem_se} from '../modules/cutadapt/se_cutadapt'
*/

include {Cutadapt} from '../modules/cutadapt/cutadapt'
include {Star_RSEM} from '../subworkflows/Star_RSEM'


workflow Mouse_RNA {

    star_genomeIndex        = Channel.of(file(params.star_genome_index, checkIfExists:true))
    rsemIndex               = Channel.of(file(params.rsem_index, checkIfExists:true))
    gtf                     = Channel.of(file(params.gtf, checkIfExists:true))
    gtf_sorted               = Channel.of(file(params.gtf_sorted, checkIfExists:true))
    gtf_sorted_index         = Channel.of(file(params.gtf_sorted_index, checkIfExists:true))


take:
    mouse_samplesheet

main:

samples_rnaseq = mouse_samplesheet
.splitCsv(header:true)
.filter { row -> row.type == "tumor_RNA" || row.type == "cell_line_RNA"}
.map { row ->
    def meta = [:]
    meta.id    =  row.sample
    meta.lib   =  row.library
    meta.sc    =  row.sample_captures
    meta.casename  = row.casename
    meta.type     = row.type
    meta.diagnosis =row.Diagnosis
    def fastq_meta = []
    fastq_meta = [ meta,  file(row.read1), file(row.read2) ]

    return fastq_meta

}

/*
Star_se(samples_rnaseq
            .combine(star_genomeIndex)
            .combine(gtf)
        )

Strandedness(
            Star_se.out.genome_bam.join(Star_se.out.genome_bai, by: [0])
            .combine(gtf_sorted)
            .combine(gtf_sorted_index)
    )

Rsem_se(Star_se.out.transcriptome_bam.join(Strandedness.out, by: [0])
            .combine(rsemIndex)
    )
*/

Cutadapt(samples_rnaseq)

Star_RSEM(Cutadapt.out.trim_reads)

}
