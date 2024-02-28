
include {Optitype
        HLA_HD
        Merge_new_HLA} from '../modules/neoantigens/hla_calls'

workflow HLA_calls {
    take: trimmed_fq
    main:

        Optitype(trimmed_fq)

        HLA_HD(trimmed_fq)

        Merge_new_HLA(Optitype.out.Optitype_output.join(HLA_HD.out.hlahd_output, by:0))

        ch_versions = Optitype.out.versions.mix(HLA_HD.out.versions)

    emit:
        mergehla = Merge_new_HLA.out
        version = ch_versions
}
