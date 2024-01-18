

include {HLAminer} from '../modules/neoantigens/hlaminer'
include {Seq2HLA} from '../modules/neoantigens/seq2hla'
include {MergeHLA} from '../modules/neoantigens/mergeHLA'
include {Optitype
        HLA_HD
        Merge_new_HLA} from '../modules/neoantigens/hla_calls'

workflow HLA_calls {
    take: trimmed_fq
    main:

        HLAminer(trimmed_fq)
        Seq2HLA(trimmed_fq)

        Optitype(trimmed_fq)

        HLA_HD(trimmed_fq)

        Merge_new_HLA(Optitype.out.Optitype_output.join(HLA_HD.out.hlahd_output, by:0))



        ch_versions = Seq2HLA.out.versions.mix(HLAminer.out.versions)
        MergeHLA(Seq2HLA.out.seq2hla_output.combine(HLAminer.out.hlaminer_output, by:0))
    emit:
        mergehla = MergeHLA.out
        version = ch_versions
}
