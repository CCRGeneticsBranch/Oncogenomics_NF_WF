include {Hotspot_Coverage} from '../modules/qc/plots'
include {Hotspot_Boxplot} from '../modules/qc/plots'
include {Flagstat} from '../modules/qc/plots'
include {Bamutil} from '../modules/qc/plots'
include {HotspotPileup} from '../modules/qc/plots'
include {Bam2tdf} from '../modules/qc/plots'
include {Coverage} from  '../modules/qc/plots'
include {CoveragePlot} from  '../modules/qc/plots'

workflow QC_from_finalBAM {

    genome                  = Channel.of(file(params.genome, checkIfExists:true))
    genome_fai              = Channel.of(file(params.genome_fai, checkIfExists:true))
    genome_dict             = Channel.of(file(params.genome_dict, checkIfExists:true))
    chrom_sizes             = Channel.of(file(params.chrom_sizes, checkIfExists:true))
    hg19_hotspot_pos         = Channel.of(file(params.hg19_hotspot_pos, checkIfExists:true))
    access_hotspot           = Channel.of(file(params.access_hotspot, checkIfExists:true))
    access_target            = Channel.of(file(params.access_target, checkIfExists:true))
    take: GATK_RNASeq_BR_PR_bam
          target_capture

    main:

     Hotspot_Coverage(
        GATK_RNASeq_BR_PR_bam
             .combine(chrom_sizes)
             .combine(access_hotspot)
     )
     Flagstat(GATK_RNASeq_BR_PR_bam)
     Bamutil(
	GATK_RNASeq_BR_PR_bam
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
     )
     HotspotPileup(
        GATK_RNASeq_BR_PR_bam
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
             .combine(hg19_hotspot_pos)
     )
     /*
     Bam2tdf(
	GATK_RNASeq_BR_PR_bam
             .combine(genome)
             .combine(genome_fai)
             .combine(genome_dict)
     )
     */
     Coverage(
        GATK_RNASeq_BR_PR_bam.combine(target_capture, by:[0])
     )
     CoveragePlot(Coverage.out.coverage_out)
     Hotspot_Boxplot(Hotspot_Coverage.out)

    emit:
         hotspot_pileup = HotspotPileup.out
         coverageplot = CoveragePlot.out
         flagstat_version = Flagstat.out.versions
         bamutil_version = Bamutil.out.versions
         flagstat = Flagstat.out.flagstat

}
