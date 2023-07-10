include {HLAminer_exome} from '../modules/neoantigens/hlaminer'
include {Seq2HLA_exome} from '../modules/neoantigens/seq2hla'
include {MergeHLA} from '../modules/neoantigens/mergeHLA'


workflow HLA_calls_exome {
    take: samples_exome
    main:
        
        HLAminer_exome(samples_exome)
        Seq2HLA_exome(samples_exome)
        MergeHLA(Seq2HLA_exome.out.combine(HLAminer_exome.out, by:0))
    emit:
        MergeHLA.out
        hlaminer_exome  = HLAminer_exome.out
        seq2hla_exome = Seq2HLA_exome.out
}
