process Picard_AddReadgroups {
        tag { dataset_id }

        input:

        tuple val(dataset_id),
        path(T_bam),
        path(G_bam)

        output:
        tuple val("${dataset_id}"),	
        path("trim_${dataset_id}_star.bam"),
        path("trim_${dataset_id}_star.bam.bai")

        container 'nciccbr/ccrgb_qctools:latest'

        script:
        """
#        java  -jar /opt2/picardtools/picardcloud.jar  AddOrReplaceReadGroups VALIDATION_STRINGENCY=SILENT INPUT=$G_bam         OUTPUT=trim_${dataset_id}_star.bam         SORT_ORDER=coordinate         RGLB=trim_${dataset_id}         RGPU=trim_${dataset_id}         RGPL=ILLUMINA         RGSM=trim_${dataset_id}	        RGCN=khanlab
       
        java  -jar /opt2/picardtools/picardcloud.jar AddOrReplaceReadGroups -VALIDATION_STRINGENCY SILENT -INPUT $G_bam  -OUTPUT trim_${dataset_id}_star.bam -SORT_ORDER coordinate -RGLB trim_${dataset_id} -RGPU trim_${dataset_id} -RGPL ILLUMINA -RGSM trim_${dataset_id} -RGCN khanlab

        java  -jar /opt2/picardtools/picardcloud.jar BuildBamIndex  -INPUT trim_${dataset_id}_star.bam -OUTPUT trim_${dataset_id}_star.bam.bai
        """
}

process Picard_MarkDuplicates {
        tag { dataset_id }

        input:

        tuple val(dataset_id),
        path(bam),
        path(index)

        output:
        tuple val("${dataset_id}"),
        path("trim_${dataset_id}.star.dd.bam"),
        path("trim_${dataset_id}.star.dd.bam.bai"),
        path("trim_${dataset_id}.markdup.txt")

        container 'nciccbr/ccrgb_qctools:latest'

        script:
        """

        java  -jar /opt2/picardtools/picardcloud.jar MarkDuplicates AS=true M=trim_${dataset_id}.markdup.txt INPUT=$bam OUTPUT=trim_${dataset_id}.star.dd.bam REMOVE_DUPLICATES=false VALIDATION_STRINGENCY=SILENT

        java  -jar /opt2/picardtools/picardcloud.jar BuildBamIndex  -INPUT trim_${dataset_id}.star.dd.bam -OUTPUT trim_${dataset_id}.star.dd.bam.bai

        """
}


