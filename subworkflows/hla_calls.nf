

include {HLAminer} from '../modules/neoantigens/hlaminer'
include {Seq2HLA} from '../modules/neoantigens/seq2hla'
include {MergeHLA} from '../modules/neoantigens/mergeHLA'


workflow HLA_calls {
    take: trimmed_fq
    main:

        HLAminer(trimmed_fq)
        Seq2HLA(trimmed_fq)
        ch_versions = Seq2HLA.out.versions.mix(HLAminer.out.versions)
        MergeHLA(Seq2HLA.out.seq2hla_output.combine(HLAminer.out.hlaminer_output, by:0))
    emit:
        mergehla = MergeHLA.out
        version = ch_versions
}
