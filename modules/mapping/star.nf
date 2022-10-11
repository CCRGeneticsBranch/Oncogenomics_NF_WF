process star {
        tag { dataset_id }
        input:
        tuple val(dataset_id),
        path(forward), 
        path(reverse),
        path(genomeIndex),
        path(gtf)
        
        output:
        tuple val("${dataset_id}"),
            path("trim_${dataset_id}Aligned.toTranscriptome.out.bam")

        container 'nciccbr/ncigb_star_v2.7.10a:latest'

        script:
        """
        STAR --genomeDir ${genomeIndex} \
                --readFilesIn $forward $reverse \
                --sjdbGTFfile ${gtf} \
                --runThreadN ${task.cpus} \
                --twopassMode Basic \
                --outSAMunmapped Within \
                --outFileNamePrefix trim_$dataset_id \
                --chimSegmentMin 12 \
                --chimJunctionOverhangMin 12 \
                --alignSJDBoverhangMin 10 \
                --alignMatesGapMax 100000 \
                --chimSegmentReadGapMax 3 \
                --outFilterMismatchNmax 2 \
                --outSAMtype BAM SortedByCoordinate \
                --quantMode TranscriptomeSAM \
                --outBAMsortingThreadN 6 \
                --limitBAMsortRAM 80000000000

        """

}

