include {Optitype
        HLA_HD
        Merge_new_HLA} from '../modules/neoantigens/hla_calls'


workflow HLA_calls_exome {
    take: samples_exome
    main:

        Optitype(samples_exome.map{ meta, R1, R2 -> [meta,[R1,R2]] })

        HLA_HD(samples_exome.map{ meta, R1, R2 -> [meta,[R1,R2]] })

        Merge_new_HLA(Optitype.out.Optitype_output.join(HLA_HD.out.hlahd_output, by:0))

        ch_versions = Optitype.out.versions.mix(HLA_HD.out.versions)

    emit:
        mergehla_exome = Merge_new_HLA.out
        ch_versions = ch_versions
}
