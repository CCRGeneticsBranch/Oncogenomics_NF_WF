process Picard_AddReadgroups {
        tag "$meta.lib"

        input:
        tuple val(meta), path(bam),path(index)

        output:
        tuple val(meta),
        path("trim_${meta.lib}_star.bam"),
        path("trim_${meta.lib}_star.bam.bai")


        script:
        def prefix = task.ext.prefix ?: "${meta.lib}"
        """
        set -exo pipefail
        printenv
          
        java -Xmx10g  -jar \$PICARDJAR AddOrReplaceReadGroups -VALIDATION_STRINGENCY SILENT -INPUT $bam  -OUTPUT trim_${prefix}_star.bam -SORT_ORDER coordinate -RGLB trim_${prefix} -RGPU trim_${prefix} -RGPL ILLUMINA -RGSM trim_${prefix} -RGCN khanlab

        java -Xmx10g  -jar \$PICARDJAR BuildBamIndex  -INPUT trim_${prefix}_star.bam -OUTPUT trim_${prefix}_star.bam.bai
        """
}

process Picard_CollectRNAseqmetrics {

        tag "$meta.lib"

       publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

        input:
        tuple val(meta),
        path(bam),
        path(index),
        path(strandedness),
        path(ref_flat),
        path(rRNA_interval)

        output:
        tuple val(meta),
        path("${meta.lib}.RnaSeqMetrics.txt"),
        path("${meta.lib}.RnaSeqMetrics.pdf")


        script:
        def prefix = task.ext.prefix ?: "${meta.lib}"
        """
        STRAND=`strandedness.py ${prefix}_strandedness.txt picard`
        java -Xmx10g -jar \$PICARDJAR CollectRnaSeqMetrics STRAND_SPECIFICITY=\$STRAND VALIDATION_STRINGENCY=SILENT REF_FLAT=$ref_flat  RIBOSOMAL_INTERVALS=$rRNA_interval INPUT=$bam OUTPUT=${prefix}.RnaSeqMetrics.txt CHART_OUTPUT=${prefix}.RnaSeqMetrics.pdf

        """
}

process Picard_CollectAlignmentSummaryMetrics {

        tag "$meta.lib"

        publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

        input:

        tuple val(meta),
        path(bam),
        path(index),
        path(genome)

        output:
        tuple val(meta),
        path("${meta.lib}.AlignmentSummaryMetrics.txt")


        script:
        """
	java  -Xmx10g -jar \$PICARDJAR CollectAlignmentSummaryMetrics VALIDATION_STRINGENCY=SILENT REFERENCE_SEQUENCE=$genome INPUT=$bam OUTPUT=${meta.lib}.AlignmentSummaryMetrics.txt ADAPTER_SEQUENCE=null

        """
}


process Picard_MarkDuplicates {
        tag "$meta.lib"

        input:

        tuple val(meta),
        path(bam),
        path(index)

        output:
        tuple val(meta),
        path("trim_${meta.lib}.dd.bam"),
        path("trim_${meta.lib}.dd.bam.bai")


        script:
        def prefix = task.ext.prefix ?: "${meta.lib}"
        """
        java  -jar -Xmx10g \$PICARDJAR MarkDuplicates AS=true M=trim_${prefix}.markdup.txt INPUT=$bam OUTPUT=trim_${prefix}.dd.bam REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT

        java  -jar -Xmx10g \$PICARDJAR BuildBamIndex  -INPUT trim_${prefix}.dd.bam -OUTPUT trim_${prefix}.dd.bam.bai

        """
}


process RNAlibrary_customQC {
        tag "$meta.lib"

        publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

        input:

        tuple val(meta),
        path(picard_rnametrics_txt),
        path(picard_rnametrics_pdf),
        path(picard_alignmentsummarymetrics),
        path(fastqc)

        output:
        tuple val(meta),
        path("${meta.lib}.RnaSeqQC.txt")


        script:
        def prefix = task.ext.prefix ?: "${meta.lib}"
        """
        rnaseqQC.pl ${fastqc}/${prefix}_R1.trim_fastqc/fastqc_data.txt ${picard_alignmentsummarymetrics} ${picard_rnametrics_txt} ${meta.id} ${prefix} ${meta.diagnosis} > ${prefix}.RnaSeqQC.txt
        
        """
}

process RNAqc_TrancriptCoverage {
        tag "$meta.lib"

        publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}/qc", mode: "${params.publishDirMode}"

        input:

        tuple val(meta),
        path(RNA_customQC),
 
        output:
        tuple val(meta),
        path("${meta.lib}.RnaSeqQC.txt")


        script:
        def prefix = task.ext.prefix ?: "${meta.lib}"
        """
        rnaseqQC.pl ${fastqc}/${prefix}_R1.trim_fastqc/fastqc_data.txt ${picard_alignmentsummarymetrics} ${picard_rnametrics_txt} ${meta.id} ${prefix} ${meta.diagnosis} > ${prefix}.RnaSeqQC.txt
        
        """
}
