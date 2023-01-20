process Hotspot_Coverage {

//   this rule uses an older version of bedtools to generate output similar to ngs_pipeline_4.2
     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/qc", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
        path(bam),
        path(index),
        path(chrom_sizes),
        path(access_hotspot)

     output:
     tuple val("${dataset_id}"),
        path("${dataset_id}.star.hotspot.depth")

     stub:
     """
     touch "${dataset_id}.star.hotspot.depth"
     """

     shell:
     '''
     set -exo pipefail

     slopBed -i !{access_hotspot}  -g !{chrom_sizes} -b 50 > !{dataset_id}_Region.bed

     ccbr_bam_filter_by_mapq.py -i !{bam} -o !{dataset_id}_filter.bam -q 30

     /opt2/bedtools2/bin/bedtools coverage -abam !{dataset_id}_filter.bam  -b !{access_hotspot} > !{dataset_id}.star.hotspot.depth

     '''
}


process Hotspot_Boxplot {

     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/qc", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
        path(hotspot)

     output:
     tuple val("${dataset_id}"),
        path("${dataset_id}.hotspot_coverage.png")

     stub:
     """
     touch "${dataset_id}.hotspot_coverage.png"
     """

     shell:
     '''
     set -exo pipefail
     boxplot.R $PWD/ !{dataset_id}.hotspot_coverage.png !{dataset_id} 

     '''
}


process Flagstat {

     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/qc", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
        path(bam),
        path(index)

     output:
     tuple val("${dataset_id}"),
        path("${dataset_id}.star.flagstat.txt")

     stub:
     """
     touch "${dataset_id}.star.flagstat.txt"
     """

     shell:
     '''
     set -exo pipefail
     samtools flagstat !{bam} > !{dataset_id}.star.flagstat.txt

     '''
}

process Bamutil {

     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/qc", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict)

     output:
     tuple val("${dataset_id}"),
        path("${dataset_id}.star.final.squeeze.bam")

     stub:
     """
     touch "${dataset_id}.star.final.squeeze.bam"
     """

     shell:
     '''
     set -exo pipefail
     bam squeeze --in !{bam} --out !{dataset_id}.star.final.squeeze.bam --refFile !{genome}  --rmTags "PG:Z;RG:Z;BI:Z;BD:Z"
     samtools index !{dataset_id}.star.final.squeeze.bam

     '''
}


process HotspotPileup {

     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/qc", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(hg19_hotspot_pos)
//        path(hotspot_mpileup_script)

     output:
     tuple val("${dataset_id}"),
        path("${dataset_id}.star.pileup.txt")

     stub:
     """
     touch "${dataset_id}.star.pileup.txt"
     """

     shell:
     '''
     set -exo pipefail
     hotspot_mpileup.pl !{hg19_hotspot_pos} !{genome} !{bam} !{dataset_id} RNAseq access > !{dataset_id}.star.pileup.txt

     '''
}


process Bam2tdf {

     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/qc", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict)

     output:
     tuple val("${dataset_id}"),
        path("${dataset_id}.star.final.bam.tdf")

     stub:
     """
     touch "${dataset_id}.star.final.bam.tdf"
     """

     shell:
     '''
     set -exo pipefail
     igvtools count !{bam} !{dataset_id}.star.final.bam.tdf  !{genome}

     '''
}

process CoveragePlot {
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/qc", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        path(bam),
        path(index),
        path(access_target)

    output:
    tuple val("${dataset_id}"),
       path("${dataset_id}.star.coverage.txt")

    stub:
     """
     touch "${dataset_id}.star.coverage.txt"
     """

    shell:
     '''
     bedtools coverage -abam !{bam} -b  !{access_target} -hist |grep "^all" > !{dataset_id}.star.coverage.txt
     coverage.R  $PWD !{dataset_id}.coveragePlot.png !{dataset_id}     
     '''

}


