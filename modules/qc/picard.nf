process Picard_AddReadgroups {
        tag { dataset_id }

//        publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/qc", mode: "${params.publishDirMode}"

        input:

        tuple val(dataset_id),
        val(library),
        path(bam),
        path(index)

        output:
        tuple val("${dataset_id}"),	
        val("${library}"),
        path("trim_${library}_star.bam"),
        path("trim_${library}_star.bam.bai")


        script:
        """
        set -exo pipefail
        printenv
          
        java -Xmx10g  -jar \$PICARDJAR AddOrReplaceReadGroups -VALIDATION_STRINGENCY SILENT -INPUT $bam  -OUTPUT trim_${library}_star.bam -SORT_ORDER coordinate -RGLB trim_${library} -RGPU trim_${library} -RGPL ILLUMINA -RGSM trim_${library} -RGCN khanlab

        java -Xmx10g  -jar \$PICARDJAR BuildBamIndex  -INPUT trim_${library}_star.bam -OUTPUT trim_${library}_star.bam.bai
        """
}

process Picard_CollectRNAseqmetrics {

        tag { dataset_id }

       publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/${library}/qc", mode: "${params.publishDirMode}"

        input:

        tuple val(dataset_id),
        val(library),
        path(bam),
        path(index),
        path(strandedness),
        path(ref_flat),
        path(rRNA_interval)

        output:
        tuple val("${dataset_id}"),
        val("${library}"),
        path("${library}.RnaSeqMetrics.txt"),
        path("${library}.RnaSeqMetrics.pdf")


        script:
        """
        STRAND=`strandedness.py ${library}_strandedness.txt picard`
        java -Xmx10g -jar \$PICARDJAR CollectRnaSeqMetrics STRAND_SPECIFICITY=\$STRAND VALIDATION_STRINGENCY=SILENT REF_FLAT=$ref_flat  RIBOSOMAL_INTERVALS=$rRNA_interval INPUT=$bam OUTPUT=${library}.RnaSeqMetrics.txt CHART_OUTPUT=${library}.RnaSeqMetrics.pdf

        """
}

process Picard_CollectAlignmentSummaryMetrics {

        tag { dataset_id }

        publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/${library}/qc", mode: "${params.publishDirMode}"

        input:

        tuple val(dataset_id),
        val(library),
        path(bam),
        path(index),
        path(genome)

        output:
        tuple val("${dataset_id}"),
        val("${library}"),
        path("${library}.AlignmentSummaryMetrics.txt")


        script:
        """

	java  -Xmx10g -jar \$PICARDJAR CollectAlignmentSummaryMetrics VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE=$genome INPUT=$bam OUTPUT=${library}.AlignmentSummaryMetrics.txt ADAPTER_SEQUENCE=null

        """
}


process Picard_MarkDuplicates {
        tag { dataset_id }

//        publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/qc", mode: "${params.publishDirMode}"

        input:

        tuple val(dataset_id),
        val(library),
        path(bam),
        path(index)

        output:
        tuple val("${dataset_id}"),
        val("${library}"),
        path("trim_${library}.star.dd.bam"),
        path("trim_${library}.star.dd.bam.bai")


        script:
        """

        java  -jar -Xmx10g \$PICARDJAR MarkDuplicates AS=true M=trim_${library}.markdup.txt INPUT=$bam OUTPUT=trim_${library}.star.dd.bam REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT

        java  -jar -Xmx10g \$PICARDJAR BuildBamIndex  -INPUT trim_${library}.star.dd.bam -OUTPUT trim_${library}.star.dd.bam.bai

        """
}


