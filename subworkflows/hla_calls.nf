

include {HLAminer} from '../modules/neoantigens/hlaminer'
include {Seq2HLA} from '../modules/neoantigens/seq2hla'
include {MergeHLA} from '../modules/neoantigens/mergeHLA'


workflow HLA_calls {
    take: trimmed_fq
    main:
        
        HLAminer(trimmed_fq)
        Seq2HLA(trimmed_fq)
        MergeHLA(Seq2HLA.out.combine(HLAminer.out, by:0))
    emit:
        MergeHLA.out
}
