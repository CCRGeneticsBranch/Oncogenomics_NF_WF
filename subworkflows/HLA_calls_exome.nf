include {HLAminer_exome} from '../modules/neoantigens/hlaminer'
include {Seq2HLA_exome} from '../modules/neoantigens/seq2hla'
include {MergeHLA} from '../modules/neoantigens/mergeHLA'


workflow HLA_calls_exome {
    take: samples_exome
    main:

        HLAminer_exome(samples_exome)

        Seq2HLA_exome(samples_exome)

        ch_versions = HLAminer_exome.out.versions.mix(Seq2HLA_exome.out.versions)

        MergeHLA(Seq2HLA_exome.out.seq2hla_output.combine(HLAminer_exome.out.hlaminer_output, by:0))
    emit:
        mergehla_exome = MergeHLA.out
        hlaminer_exome  = HLAminer_exome.out.hlaminer_output
        seq2hla_exome = Seq2HLA_exome.out.seq2hla_output
        ch_versions = ch_versions
}
