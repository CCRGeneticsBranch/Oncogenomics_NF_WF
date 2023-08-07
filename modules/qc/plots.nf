process Hotspot_Coverage {

//   this rule uses an older version of bedtools to generate output similar to ngs_pipeline_4.2
     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
        path(bam),
        path(index),
        path(chrom_sizes),
        path(access_hotspot)

     output:
     tuple val(meta),
        path("${meta.lib}.hotspot.depth")

     stub:
     """
     touch "${meta.lib}.hotspot.depth"
     """

    script:
      def prefix = task.ext.prefix ?: "${meta.lib}"
      """
     set -exo pipefail

     slopBed -i ${access_hotspot}  -g ${chrom_sizes} -b 50 > ${prefix}_Region.bed

     ccbr_bam_filter_by_mapq.py -i ${bam} -o ${prefix}_filter.bam -q 30

     /opt2/bedtools2/bin/bedtools coverage -abam ${prefix}_filter.bam  -b ${access_hotspot} > ${prefix}.hotspot.depth

     """
}


process Hotspot_Boxplot {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/qc", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
        path(hotspot)

     output:
     tuple val(meta),
        path("${meta.id}.hotspot_coverage.png")

     stub:
     """
     touch "${meta.id}.hotspot_coverage.png"
     """

     script:
     """
     set -exo pipefail
     boxplot.R \$PWD/ ${meta.id}.hotspot_coverage.png ${meta.lib} 

     """
}

process Flagstat {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
        path(bam),
        path(index)

     output:
     tuple val(meta),
        path("${meta.lib}.flagstat.txt")

     stub:
     """
     touch "${meta.lib}.flagstat.txt"
     """

     script:
     """
     set -exo pipefail
     samtools flagstat ${bam} > ${meta.lib}.flagstat.txt

     """
}


process Bamutil {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict)

     output:
     tuple val(meta),
        path("${meta.lib}.final.squeeze.bam")

     stub:
     """
     touch "${meta.lib}.final.squeeze.bam"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
      """
     set -exo pipefail
     bam squeeze --in ${bam} --out ${prefix}.final.squeeze.bam --refFile ${genome}  --rmTags "PG:Z;RG:Z;BI:Z;BD:Z"
     samtools index ${prefix}.final.squeeze.bam
    
     """
}


process HotspotPileup {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(hg19_hotspot_pos)


     output:
     tuple val(meta),
        path("${meta.lib}.pileup.txt")

     stub:
     """
     touch "${meta.lib}.pileup.txt"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
     set -exo pipefail
     hotspot_mpileup.pl ${hg19_hotspot_pos} ${genome} ${bam} ${prefix} ${meta.type} ${meta.sc} > ${prefix}.pileup.txt

     """
}

process MakeHotSpotDB {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.id}/db", mode: "${params.publishDirMode}"

     input:
     path(libs)
     val(meta)
     //tuple val(meta), path(pileup)
   
     output:
     tuple val(meta),
        path("${meta.id}.hotspot")

     stub:
     """
     touch "${meta.id}.hotspot"
     """

     script:
 
     """
     cat ${libs.join(' ')} |sort > ${meta.id}.hotspot
     """
}


process Bam2tdf {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}", mode: "${params.publishDirMode}"

     input:
     tuple val(meta),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict)

     output:
     tuple val(meta),
        path("${meta.lib}.final.bam.tdf")

     stub:
     """
     touch "${meta.lib}.final.bam.tdf"
     """

     script:
     """
     set -exo pipefail
     igvtools count ${bam} ${meta.lib}.final.bam.tdf  ${genome}

     """
}

process Coverage {
   
    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(bam),
        path(index),
        path(targetcapture)

    output:
    tuple val(meta),
       path("${meta.lib}.coverage.txt")

    stub:
     """
     touch "${meta.lib}.coverage.txt"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
     bedtools coverage -abam ${bam} -b  ${targetcapture} -hist |grep "^all" > ${prefix}.coverage.txt
         
     """

}

process CoveragePlot {
   
    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(coverage)

    output:
    tuple val(meta),
       path("${meta.lib}.coveragePlot.png")

    stub:
     """
     touch "${meta.lib}.coveragePlot.png"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
     coverage.R  \$PWD ${prefix}.coveragePlot.png ${prefix}

     """

}

process Read_depth {
   
    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(bam),
        path(index),
        path(targetcapture)

   output:
    tuple val(meta),
       path("${meta.lib}.depth_per_base")
   stub:
     """
     touch "${meta.lib}.depth_per_base"
     """
   script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
   """
   echo -e "chr\tstart\tend\tgene\tposition\tdepth" >  ${prefix}.depth_per_base
   cut -f1-4 ${targetcapture} > intervals.bed
   samtools view -hF 0x400 -q 30 -L intervals.bed ${bam} |samtools view -ShF 0x4 - | samtools view -SuF 0x200 - | bedtools coverage -split -a intervals.bed -b - -d >> ${prefix}.depth_per_base

   """

}

process VerifyBamID {
    
    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/verifyBamID", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(bam),
        path(index),
        path(recode_vcf)

    output:
    tuple val(meta),
        path("${meta.lib}.selfSM")

    stub:
    """
    touch "${meta.lib}.selfSM"
    """

    script:

     """
     verifyBamID --vcf ${recode_vcf} --bam ${bam} --maxDepth 3000 --ignoreRG --site --chip-none --precise --minMapQ 30 --minQ 20 --minAF 0.05 --out \$PWD/${meta.lib}
     """
}

process FailedExons_Genes {
   
    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(depth_per_base)

   output:
    tuple val(meta),
       path("${meta.lib}.failExons"),
       path("${meta.lib}.failGenes")

   stub:
     """
     touch "${meta.lib}.failExons"
     touch "${meta.lib}.failGenes"
     """
   script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
   """
   failed_Exon_Final.pl ${depth_per_base} 10 ${prefix}.failExons ${prefix}.failGenes
   """
}


process TargetIntervals {
   
    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(bam),
        path(index),
        path(targetcapture),
        path(design_ch)

    output:
    tuple val(meta),
       path("${meta.lib}.probe.intervals"),
       path("${meta.lib}.target.intervals")

    stub:
     """
     touch "${meta.lib}.probe.intervals"
     touch "${meta.lib}.target.intervals"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
   
     cat <(samtools view -H ${bam}) <(awk '{{print \$1 "\t" \$2+1 "\t" \$3 "\t+\tinterval_" NR}}' ${design_ch} )> ${prefix}.probe.intervals
     cat <(samtools view -H ${bam}) <(awk '{{print \$1 "\t" \$2+1 "\t" \$3 "\t+\tinterval_" NR}}' ${targetcapture} )> ${prefix}.target.intervals
     """

}

process HSMetrics {

    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(bam),
        path(index),
        path(probe_intervals),
        path(target_intervals),
        path(genome),
        path(genome_fai),
        path(genome_dict)

    output:
    tuple val(meta),
       path("${meta.lib}.hsmetrics")

    stub:
     """
     touch "${meta.lib}.hsmetrics"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
     java -Xmx75g  -jar \$PICARDJAR CollectHsMetrics BAIT_INTERVALS= ${probe_intervals} TARGET_INTERVALS= ${target_intervals} INPUT= ${bam} OUTPUT= ${prefix}.hsmetrics METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE= ${genome} QUIET=true  VALIDATION_STRINGENCY=SILENT
     """

}

