process Star {
    tag { dataset_id }

    publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/${library}", mode: "${params.publishDirMode}"

    input:
    tuple val(dataset_id),
        val(library),
        path(r1), 
        path(r2),
        path(star_genomeIndex),
        path(gtf)
    
    output:
    tuple val("${dataset_id}"), val("${library}"), path("${library}.Aligned.toTranscriptome.out.bam") , emit: transcriptome_bam
    tuple val("${dataset_id}"), val("${library}"), path("${library}.Aligned.sortedByCoord.out.bam"), emit: genome_bam
    tuple val("${dataset_id}"), val("${library}"), path("${library}.Aligned.sortedByCoord.out.bam.bai"), emit: genome_bai
    tuple val("${dataset_id}"), val("${library}"), path("${library}.Chimeric.out.junction"), emit: chimeric_junction

    stub:
    """
    touch "${library}.Aligned.toTranscriptome.out.bam"
    touch "${library}.Aligned.sortedByCoord.out.bam"
    touch "${library}.Aligned.sortedByCoord.out.bam.bai"
    touch "${library}.Chimeric.out.junction"
    """

    shell:
    '''
    set -exo pipefail
    if [ -d "/lscratch/${SLURM_JOB_ID}" ];then
        TMPDIR="/lscratch/${SLURM_JOB_ID}/!{library}_STAR"
        if [ -d ${TMPDIR} ];then rm -rf ${TMPDIR};fi
        
        # run STAR alignment
        STAR --genomeDir !{star_genomeIndex} \
            --readFilesIn !{r1} !{r2} \
            --readFilesCommand zcat \
            --sjdbGTFfile !{gtf} \
            --runThreadN !{task.cpus} \
            --twopassMode Basic \
            --outSAMunmapped Within \
            --outFileNamePrefix "!{library}." \
            --chimSegmentMin 12 \
            --chimOutJunctionFormat 1 \
            --chimJunctionOverhangMin 12 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 100000 \
            --chimSegmentReadGapMax 3 \
            --outFilterMismatchNmax 2 \
            --outSAMtype BAM Unsorted \
            --outTmpDir ${TMPDIR} \
            --quantMode TranscriptomeSAM 
            
        # sort files
        samtools sort -@ !{task.cpus} -T ${TMPDIR} -o !{library}.Aligned.sortedByCoord.out.bam -O BAM !{library}.Aligned.out.bam
    
    else
        # run STAR alignment
        STAR --genomeDir !{star_genomeIndex} \
            --readFilesIn !{r1} !{r2} \
            --readFilesCommand zcat \
            --sjdbGTFfile !{gtf} \
            --runThreadN !{task.cpus} \
            --twopassMode Basic \
            --outSAMunmapped Within \
            --outFileNamePrefix "!{library}." \
            --chimSegmentMin 12 \
            --chimOutJunctionFormat 1 \
            --chimJunctionOverhangMin 12 \
            --alignSJDBoverhangMin 10 \
            --alignMatesGapMax 100000 \
            --chimSegmentReadGapMax 3 \
            --outFilterMismatchNmax 2 \
            --outSAMtype BAM Unsorted \
            --quantMode TranscriptomeSAM 
            
        # sort files
        samtools sort -@ !{task.cpus} -o !{library}.Aligned.sortedByCoord.out.bam -O BAM !{library}.Aligned.out.bam
    fi
    
    # index files
    samtools index -@ !{task.cpus} !{library}.Aligned.sortedByCoord.out.bam
    '''
}
