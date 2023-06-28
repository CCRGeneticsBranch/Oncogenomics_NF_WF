include {Exome_common_WF} from './Exome_common_WF.nf'



workflow Tumor_Normal_WF {
  

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
    def fastq_meta = []
    fastq_meta = [ meta,  file(row.read1), file(row.read2)  ]

    return fastq_meta
}
//samples_exome.view()


Exome_common_WF(samples_exome)


}