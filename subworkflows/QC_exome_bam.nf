
include {Read_depth} from  '../modules/qc/plots'
include {Bam2tdf} from '../modules/qc/plots'
include {Flagstat} from '../modules/qc/plots'
include {Bamutil} from '../modules/qc/plots'
include {HotspotPileup} from '../modules/qc/plots'
include {Genotyping} from '../modules/qc/qc'
include {CircosPlot_lib} from '../modules/qc/qc'
include {VerifyBamID} from '../modules/qc/plots'
include {FailedExons_Genes} from '../modules/qc/plots'
include {Coverage} from  '../modules/qc/plots'
include {CoveragePlot} from  '../modules/qc/plots'
include {TargetIntervals} from  '../modules/qc/plots'
include {HSMetrics} from  '../modules/qc/plots'
include {Conpair_pile} from  '../modules/qc/qc'
include {Exome_QC} from '../modules/qc/qc.nf'
include {Hotspot_Coverage} from '../modules/qc/plots'

workflow QC_exome_bam {

    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
    chrom_sizes             = Channel.of(file(params.chrom_sizes, checkIfExists:true))
    hg19_hotspot_pos        = Channel.of(file(params.hg19_hotspot_pos, checkIfExists:true))
    clin_ex_v1              = Channel.of(file(params.clin_ex_v1, checkIfExists:true))
    Sites1000g4genotyping   = Channel.of(file(params.Sites1000g4genotyping, checkIfExists:true))
    recode_vcf              = Channel.of(file(params.recode_vcf, checkIfExists:true))
    conpair_refbed          = Channel.of(file(params.conpair_refbed, checkIfExists:true))
    access_hotspot           = Channel.of(file(params.access_hotspot, checkIfExists:true))

    take:
         GATK_exome_bam
         bwa_picard_bam
         capture_ch
         design_ch
    main:
         Genotyping(
          bwa_picard_bam
            .combine(Sites1000g4genotyping)
            .combine(genome)
            .combine(genome_fai)
            .combine(genome_dict)
         )
         CircosPlot_lib(Genotyping.out.loh)
         Read_depth(GATK_exome_bam.combine(capture_ch, by:[0]))

         ch_versions = Genotyping.out.versions.mix(Read_depth.out.versions)

         FailedExons_Genes(Read_depth.out.read_depth_output)
         /*
         Bam2tdf(
          GATK_exome_bam
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
         )
        */
         //ch_versions = ch_versions.mix(Bam2tdf.out.versions)

         Flagstat(GATK_exome_bam)

         ch_versions = ch_versions.mix(Flagstat.out.versions)

         Bamutil(
          GATK_exome_bam
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
         )

         ch_versions = ch_versions.mix(Bamutil.out.versions)

        HotspotPileup(
          GATK_exome_bam
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
             .combine(hg19_hotspot_pos)
        )
        VerifyBamID(
          GATK_exome_bam
             .combine(recode_vcf)
        )

        ch_versions = ch_versions.mix(VerifyBamID.out.versions)

        Coverage(
        GATK_exome_bam.combine(capture_ch, by:[0])
        )

        ch_versions = ch_versions.mix(Coverage.out.versions)

        CoveragePlot(Coverage.out.coverage_out)
        TargetIntervals(
          GATK_exome_bam.combine(capture_ch, by:[0]).combine(design_ch, by:[0])
        )

        ch_versions = ch_versions.mix(TargetIntervals.out.versions)

        Conpair_pile(
          GATK_exome_bam
          .combine(genome)
          .combine(genome_fai)
          .combine(genome_dict)
          .combine(conpair_refbed)
        )

        HSMetrics(GATK_exome_bam
          .combine(TargetIntervals.out.probe_intervals, by:[0])
          .combine(TargetIntervals.out.target_intervals, by:[0])
          .combine(genome)
          .combine(genome_fai)
          .combine(genome_dict)
          )

        ch_versions = ch_versions.mix(HSMetrics.out.versions)

        Exome_QC(
          GATK_exome_bam.combine(capture_ch, by:[0]).combine(HSMetrics.out.hsmetrics_out, by:[0])
        )
        Hotspot_Coverage(
        GATK_exome_bam
             .combine(chrom_sizes)
             .combine(access_hotspot)
        )

   emit:
         hotspot_pileup = HotspotPileup.out
         coverageplot = CoveragePlot.out
         loh = Genotyping.out.loh
         gt = Genotyping.out.gt
         Exome_qc = Exome_QC.out
         hsmetrics = HSMetrics.out.hsmetrics_out
         flagstat = Flagstat.out.flagstat
         verifybamid = VerifyBamID.out.verifybamid_out
         ch_versions = ch_versions
         conpair_pileup = Conpair_pile.out

}
