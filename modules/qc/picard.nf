
process Picard_AddReadgroups {
        tag { dataset_id }

        publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/qc", mode: "${params.publishDirMode}"

        input:

        tuple val(dataset_id),
        path(bam),
        path(index)

        output:
        tuple val("${dataset_id}"),	
        path("trim_${dataset_id}_star.bam"),
        path("trim_${dataset_id}_star.bam.bai")


        script:
        """
        set -exo pipefail
        printenv
          
        java  -jar \$PICARDJAR AddOrReplaceReadGroups -VALIDATION_STRINGENCY SILENT -INPUT $bam  -OUTPUT trim_${dataset_id}_star.bam -SORT_ORDER coordinate -RGLB trim_${dataset_id} -RGPU trim_${dataset_id} -RGPL ILLUMINA -RGSM trim_${dataset_id} -RGCN khanlab

        java  -jar \$PICARDJAR BuildBamIndex  -INPUT trim_${dataset_id}_star.bam -OUTPUT trim_${dataset_id}_star.bam.bai
        """
}

process Picard_CollectRNAseqmetrics {

        tag { dataset_id }

       publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/qc", mode: "${params.publishDirMode}"

        input:

        tuple val(dataset_id),
        path(bam),
        path(index),
        path(ref_flat),
        path(rRNA_interval)

        output:
        tuple val("${dataset_id}"),
        path("trim_${dataset_id}.RnaSeqMetrics.txt"),
        path("trim_${dataset_id}.RnaSeqMetrics.pdf")


        script:
        """

        java  -jar \$PICARDJAR CollectRnaSeqMetrics STRAND_SPECIFICITY=NONE VALIDATION_STRINGENCY=SILENT REF_FLAT=$ref_flat  RIBOSOMAL_INTERVALS=$rRNA_interval INPUT=$bam OUTPUT=trim_${dataset_id}.RnaSeqMetrics.txt CHART_OUTPUT=trim_${dataset_id}.RnaSeqMetrics.pdf

        """
}

process Picard_CollectAlignmentSummaryMetrics {

        tag { dataset_id }

        publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/qc", mode: "${params.publishDirMode}"

        input:

        tuple val(dataset_id),
        path(bam),
        path(index),
        path(genome)

        output:
        tuple val("${dataset_id}"),
        path("trim_${dataset_id}.AlignmentSummaryMetrics.txt")


        script:
        """

	java  -jar \$PICARDJAR CollectAlignmentSummaryMetrics VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE=$genome INPUT=$bam OUTPUT=trim_${dataset_id}.AlignmentSummaryMetrics.txt ADAPTER_SEQUENCE=null

        """
}


process Picard_MarkDuplicates {
        tag { dataset_id }

        publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/qc", mode: "${params.publishDirMode}"

        input:

        tuple val(dataset_id),
        path(bam),
        path(index)

        output:
        tuple val("${dataset_id}"),
        path("trim_${dataset_id}.star.dd.bam"),
        path("trim_${dataset_id}.star.dd.bam.bai")


        script:
        """

        java  -jar \$PICARDJAR MarkDuplicates AS=true M=trim_${dataset_id}.markdup.txt INPUT=$bam OUTPUT=trim_${dataset_id}.star.dd.bam REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT

        java  -jar \$PICARDJAR BuildBamIndex  -INPUT trim_${dataset_id}.star.dd.bam -OUTPUT trim_${dataset_id}.star.dd.bam.bai

        """
}


