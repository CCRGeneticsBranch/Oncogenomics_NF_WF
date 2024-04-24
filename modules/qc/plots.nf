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
        path("${meta.lib}.*.hotspot.depth")

     stub:
     """
     touch "${meta.lib}.*.hotspot.depth"
     """

    script:
      def prefix = task.ext.prefix ?: "${meta.lib}"
      """
     set -exo pipefail

     slopBed -i ${access_hotspot}  -g ${chrom_sizes} -b 50 > ${prefix}_Region.bed

     ccbr_bam_filter_by_mapq.py -i ${bam} -o ${prefix}_filter.bam -q 30

     /opt2/bedtools2/bin/bedtools coverage -abam ${prefix}_filter.bam  -b ${access_hotspot} > ${prefix}.hotspot.depth

    if [[ "${meta.type}" == *"DNA"* ]]; then
        mv ${prefix}.hotspot.depth ${prefix}.bwa.hotspot.depth
    else
        mv ${prefix}.hotspot.depth ${prefix}.star.hotspot.depth
    fi

     """
}


process Hotspot_Boxplot {

     tag "$meta.id"

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
     boxplot.R \$PWD/ ${meta.id}.hotspot_coverage.png ${meta.id}

     """
}

process Flagstat {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

     input:
     tuple val(meta),
        path(bam),
        path(index)

     output:
     tuple val(meta),path("${meta.lib}.flagstat.txt"), emit: flagstat
     path "versions.yml"             , emit: versions


     stub:
     """
     touch "${meta.lib}.flagstat.txt"
     """

     script:
     """
     echo ${workflow.homeDir}
     set -exo pipefail
     samtools flagstat ${bam} > ${meta.lib}.flagstat.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version|head -n 1 |sed 's/.*samtools //')
    END_VERSIONS

     """
}


process Bamutil {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

     input:
     tuple val(meta),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict)

     output:
     tuple val(meta),
        path("${meta.lib}.final.squeeze.bam"),
        path("${meta.lib}.final.squeeze.bam.bai")
     path "versions.yml"             , emit: versions

     stub:
     """
     touch "${meta.lib}.final.squeeze.bam"
     touch "${meta.lib}.final.squeeze.bam.bai"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
      """
     set -exo pipefail
     bam squeeze --in ${bam} --out ${prefix}.final.squeeze.bam --refFile ${genome}  --rmTags "PG:Z;RG:Z;BI:Z;BD:Z"
     samtools index ${prefix}.final.squeeze.bam
     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
         samtools: \$(samtools --version|head -n 1 |sed 's/.*samtools //')
     END_VERSIONS

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

process MakeHotSpotDB_old {

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

process MakeHotSpotDB {

     tag "$meta.id"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.id}/db", mode: "${params.publishDirMode}"

     input:

    tuple val(meta), path(pile_files)

     output:
     tuple val(meta),
        path("${meta.id}.hotspot")

     stub:
     """
     touch "${meta.id}.hotspot"
     """

     script:

     """
     cat ${pile_files.join(' ')} |sort > ${meta.id}.hotspot
     """
}


process Bam2tdf {

     tag "$meta.lib"

     publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

     input:
     tuple val(meta),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict)

     output:
     tuple val(meta),path("${meta.lib}.final.bam.tdf"), emit: tdf_out
     path "versions.yml"             , emit: versions

     stub:
     """
     touch "${meta.lib}.final.bam.tdf"
     """

     script:
     """
     set -exo pipefail
     igvtools count ${bam} ${meta.lib}.final.bam.tdf  ${genome}

     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
         igvtools: \$(igvtools 2>&1|grep -E '^Program'|cut -f5 -d " ")
     END_VERSIONS
     """
}

process Coverage {

    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

    input:
    tuple val(meta),
        path(bam),
        path(index),
        path(targetcapture),
        path(sorted_chr_order),
        path(genomelength)

    output:
    tuple val(meta),path("${meta.lib}.*.coverage.txt"), emit:coverage_out
    path "versions.yml"             , emit: versions


    stub:
     """
     touch "${meta.lib}.*.coverage.txt"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
     if ! grep -q "^chrM" ${targetcapture} ; then echo -e "chrM\t3306\t15887" >> ${targetcapture} ; fi
     awk -F'\t'  '\$1 !~ /_/' ${targetcapture}|awk 'NR==FNR{order[\$1]=NR; next} {print order[\$1]"\t"\$0}' ${sorted_chr_order} - | \
     sort -k1,1n -k3,3n |tr -s ' ' '\t' | cut -f 2,3,4 |sed 's/\t*\$//' > sorted_bed
     bedtools coverage -a sorted_bed -sorted -b ${bam} -g ${genomelength} -hist |grep "^all" > ${prefix}.coverage.txt

     if [[ "${meta.type}" == *"DNA"* ]]; then
        mv ${prefix}.coverage.txt ${prefix}.bwa.coverage.txt
     else
        mv ${prefix}.coverage.txt ${prefix}.star.coverage.txt
     fi

     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
         bedtools: \$(bedtools --version|sed 's/.*bedtools //')
     END_VERSIONS

     """

}

process CoveragePlot {

    tag "$meta.id"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(meta),
        path(coverage)

    output:
    tuple val(meta),
       path("${meta.id}.coveragePlot.png")

    stub:
     """
     touch "${meta.id}.coveragePlot.png"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.id}"
     """
     TMP=tmp/
     mkdir \$TMP
     trap 'rm -rf "\$TMP"' EXIT
     cp ${coverage.join(' ')} \$TMP
     coverage.R  \$TMP ${prefix}.coveragePlot.png ${prefix}

     """

}

process Read_depth {

    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

    input:
    tuple val(meta),
        path(bam),
        path(index),
        path(targetcapture),
        path(sorted_chr_order),
        path(genomelength)

   output:
    tuple val(meta),path("${meta.lib}.depth_per_base"), emit: read_depth_output
    path "versions.yml"             , emit: versions
   stub:
     """
     touch "${meta.lib}.depth_per_base"
     """
   script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
   """
   echo -e "chr\tstart\tend\tgene\tposition\tdepth" >  ${prefix}.depth_per_base

   awk -F'\t'  '\$1 !~ /_/' ${targetcapture}|sed  '/^chrM/d'|awk 'NR==FNR{order[\$1]=NR; next} {print order[\$1]"\t"\$0}' ${sorted_chr_order} - | \
   sort -k1,1n -k3,3n |tr -s ' ' '\t' |cut -f2-  > sorted_bed
   cut -f1-4 sorted_bed > intervals.bed

   samtools view -hF 0x400 -q 30 -L intervals.bed ${bam} |samtools view -ShF 0x4 - | samtools view -SuF 0x200 - | bedtools coverage -split -a intervals.bed -sorted -b - -g ${genomelength} -d >> ${prefix}.depth_per_base

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       samtools: \$(samtools --version|grep -E '^samtools'|sed 's/.*samtools //')
   END_VERSIONS
   """

}

process VerifyBamID {

    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/verifyBamID", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

    input:
    tuple val(meta),
        path(bam),
        path(index),
        path(recode_vcf)

    output:
    tuple val(meta),path("${meta.lib}.selfSM"), emit: verifybamid_out
    path "versions.yml"             , emit: versions


    stub:
    """
    touch "${meta.lib}.selfSM"
    """

    script:

     """
     verifyBamID --vcf ${recode_vcf} --bam ${bam} --maxDepth 3000 --ignoreRG --site --chip-none --precise --minMapQ 30 --minQ 20 --minAF 0.05 --out \$PWD/${meta.lib}

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       verifyBamID: \$(verifyBamID 2>&1|grep -E '^verifyBamID'|cut -f2 -d " ")
   END_VERSIONS

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
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

    input:
    tuple val(meta),
        path(bam),
        path(index),
        path(targetcapture),
        path(design_ch)

    output:
    tuple val(meta),path("${meta.lib}.probe.intervals"), emit: probe_intervals
    tuple val(meta),path("${meta.lib}.target.intervals"), emit: target_intervals
    path "versions.yml"             , emit: versions

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

   cat <<-END_VERSIONS > versions.yml
   "${task.process}":
       samtools: \$(samtools --version|grep -E '^samtools'|sed 's/.*samtools //')
   END_VERSIONS
     """

}

process HSMetrics {

    tag "$meta.lib"
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}",pattern: "${meta.lib}*"

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
    tuple val(meta),path("${meta.lib}.hsmetrics"), emit: hsmetrics_out
    path "versions.yml"             , emit: versions

    stub:
     """
     touch "${meta.lib}.hsmetrics"
     """

     script:
     def prefix = task.ext.prefix ?: "${meta.lib}"
     """
     java -Xmx75g  -jar \$PICARDJAR CollectHsMetrics BAIT_INTERVALS= ${probe_intervals} TARGET_INTERVALS= ${target_intervals} INPUT= ${bam} OUTPUT= ${prefix}.histogram.hsmetrics METRIC_ACCUMULATION_LEVEL=ALL_READS REFERENCE_SEQUENCE= ${genome} QUIET=true  VALIDATION_STRINGENCY=SILENT
     awk '/^## HISTOGRAM/ {{exit}} {{print}}' ${prefix}.histogram.hsmetrics > ${prefix}.hsmetrics

     cat <<-END_VERSIONS > versions.yml
     "${task.process}":
         Picard: \$(java -jar \$PICARDJAR MarkDuplicates --version 2>&1 |sed 's/Version://')
     END_VERSIONS
     """

}
