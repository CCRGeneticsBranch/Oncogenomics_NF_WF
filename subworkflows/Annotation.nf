
include {Annovar} from '../modules/annotation/annot'
include {Custom_annotation} from '../modules/annotation/annot'
include {Combine_annotation} from '../modules/annotation/annot'

workflow Annotation {

    annovar_data             = Channel.of(file(params.annovar_data, checkIfExists:true))
    clinseq                  = Channel.of(file(params.clinseq, checkIfExists:true))
    cosmic                   = Channel.of(file(params.cosmic, checkIfExists:true))
    pcg                      = Channel.of(file(params.pcg, checkIfExists:true))
    clinvar                  = Channel.of(file(params.clinvar, checkIfExists:true))
    hgmd                     = Channel.of(file(params.hgmd, checkIfExists:true))
    matchTrial               = Channel.of(file(params.matchTrial, checkIfExists:true))
    mcg                      = Channel.of(file(params.mcg, checkIfExists:true))
    DoCM                     = Channel.of(file(params.DoCM, checkIfExists:true))
    CanDL                    = Channel.of(file(params.CanDL, checkIfExists:true))
    targetted_cancer_care    = Channel.of(file(params.targetted_cancer_care, checkIfExists:true))
    civic                    = Channel.of(file(params.civic, checkIfExists:true))
    ACMG                     = Channel.of(file(params.ACMG, checkIfExists:true))
    hg19_BLsites             = Channel.of(file(params.hg19_BLsites, checkIfExists:true))
    hg19_WLsites             = Channel.of(file(params.hg19_WLsites, checkIfExists:true))

    take:
        FormatInput_out
    main:

    Annovar(
        FormatInput_out
             .combine(annovar_data)
             .combine(clinseq)
             .combine(pcg)
    )
    Custom_annotation(
         FormatInput_out
             .combine(clinvar)
             .combine(hgmd)
             .combine(matchTrial)
             .combine(mcg)
             .combine(DoCM)
             .combine(CanDL)
             .combine(targetted_cancer_care)
             .combine(civic)
    )
    Annovar_out_ch = Annovar.out.cosmic
                        .join(Annovar.out.clinseq, by:[0])
                        .join(Annovar.out.cadd, by:[0])
                        .join(Annovar.out.pcg, by:[0])
                        .join(Annovar.out.gene, by:[0])
    Combine_annotation(Annovar_out_ch.join(Custom_annotation.out, by:[0])
             .combine(ACMG)
             .combine(hg19_BLsites)
             .combine(hg19_WLsites)
    )

    emit:
    rare_annotation	= Combine_annotation.out.rare_annotation
    final_annotation	= Combine_annotation.out.final_annotation
    version = Annovar.out.versions


}
