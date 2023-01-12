process Hotspot_Coverage {

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
     samtools view -hF 0x400 -q 30 -L !{dataset_id}_Region.bed !{bam} | samtools view -ShF 0x4 - | samtools view -SuF 0x200 - | bedtools coverage -abam - -b !{access_hotspot} > !{dataset_id}.star.hotspot.depth
     '''
}


process Hotspot_Boxplot {

     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/qc", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
        path(hotspot),
        path(boxplot_script)

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
     R --vanilla --slave --silent --args !{hotspot} !{dataset_id}.hotspot_coverage.png !{dataset_id} <!{boxplot_script}

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



