process Star {
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/STAR", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        path(r1), 
        path(r2),
        path(star_genomeIndex),
        path(gtf)
    
    output:
    tuple val("${dataset_id}"),
        path("${dataset_id}.Aligned.toTranscriptome.out.bam"),
        path("${dataset_id}.Aligned.sortedByCoord.out.bam"),
        path("${dataset_id}.Aligned.sortedByCoord.out.bam.bai"),
        path("${dataset_id}.Chimeric.out.junction")

    stub:
    """
    touch "${dataset_id}.Aligned.toTranscriptome.out.bam"
    touch "${dataset_id}.Aligned.sortedByCoord.out.bam"
    touch "${dataset_id}.Aligned.sortedByCoord.out.bam.bai"
    touch "${dataset_id}.Chimeric.out.junction"
    """


    shell:
    '''
    set -exo pipefail
    if [ -d "/lscratch/${SLURM_JOB_ID}" ];then
        TMPDIR="/lscratch/${SLURM_JOB_ID}/!{dataset_id}_STAR"
        if [ -d ${TMPDIR} ];then rm -rf ${TMPDIR};fi
        
        # run STAR alignment
        STAR --genomeDir !{star_genomeIndex} \
            --readFilesIn !{r1} !{r2} \
            --readFilesCommand zcat \
            --sjdbGTFfile !{gtf} \
            --runThreadN !{task.cpus} \
            --twopassMode Basic \
            --outSAMunmapped Within \
            --outFileNamePrefix "!{dataset_id}." \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 12 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 100000 \
            --chimSegmentReadGapMax 3 \
            --outFilterMismatchNmax 2 \
            --outSAMtype BAM Unsorted \
            --outTmpDir ${TMPDIR} \
            --quantMode TranscriptomeSAM 
            
        # sort files
        samtools sort -@ !{task.cpus} -T ${TMPDIR} -o !{dataset_id}.Aligned.sortedByCoord.out.bam -O BAM !{dataset_id}.Aligned.out.bam
    
    else
        # run STAR alignment
        STAR --genomeDir !{star_genomeIndex} \
            --readFilesIn !{r1} !{r2} \
            --readFilesCommand zcat \
            --sjdbGTFfile !{gtf} \
            --runThreadN !{task.cpus} \
            --twopassMode Basic \
            --outSAMunmapped Within \
            --outFileNamePrefix "!{dataset_id}." \
            --chimSegmentMin 12 \
            --chimJunctionOverhangMin 12 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 100000 \
            --chimSegmentReadGapMax 3 \
            --outFilterMismatchNmax 2 \
            --outSAMtype BAM Unsorted \
            --quantMode TranscriptomeSAM 
            
        # sort files
        samtools sort -@ !{task.cpus} -o !{dataset_id}.Aligned.sortedByCoord.out.bam -O BAM !{dataset_id}.Aligned.out.bam
    fi
    
    # index files
    samtools index -@ !{task.cpus} !{dataset_id}.Aligned.sortedByCoord.out.bam
    '''
}
