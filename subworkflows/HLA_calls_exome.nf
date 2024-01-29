include {HLAminer_exome} from '../modules/neoantigens/hlaminer'
include {Seq2HLA_exome} from '../modules/neoantigens/seq2hla'
include {MergeHLA} from '../modules/neoantigens/mergeHLA'
include {Optitype
        HLA_HD
        Merge_new_HLA} from '../modules/neoantigens/hla_calls'


workflow HLA_calls_exome {
    take: samples_exome
    main:

        //HLAminer_exome(samples_exome)

        //Seq2HLA_exome(samples_exome)

        Optitype(samples_exome.map{ meta, R1, R2 -> [meta,[R1,R2]] })

        HLA_HD(samples_exome.map{ meta, R1, R2 -> [meta,[R1,R2]] })

        Merge_new_HLA(Optitype.out.Optitype_output.join(HLA_HD.out.hlahd_output, by:0))

        ch_versions = Optitype.out.versions.mix(HLA_HD.out.versions)

        //MergeHLA(Seq2HLA_exome.out.seq2hla_output.combine(HLAminer_exome.out.hlaminer_output, by:0))
    emit:
        mergehla_exome = Merge_new_HLA.out
        //hlaminer_exome  = HLAminer_exome.out.hlaminer_output
        //seq2hla_exome = Seq2HLA_exome.out.seq2hla_output
        ch_versions = ch_versions
}
