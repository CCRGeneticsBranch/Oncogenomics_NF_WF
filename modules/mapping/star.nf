process Star {
    tag { dataset_id }

    publishDir "$params.resultsdir/$dataset_id/STAR", mode: 'copy'

    input:
    tuple val(dataset_id),
        path(r1), 
        path(r2),
        path(genomeIndex),
        path(gtf)
    
    output:
    tuple val("${dataset_id}"),
        path("${dataset_id}.Aligned.toTranscriptome.out.bam"),
        path("${dataset_id}.Aligned.sortedByCoord.out.bam"),
        path("${dataset_id}.Aligned.sortedByCoord.out.bam.bai")

    script:
    """
    set -exo pipefail
    if [ -d /lscratch/\${SLURM_JOB_ID} ];then
        TMPDIR="/lscratch/${SLURM_JOB_ID}/$dataset_id"
    else
        TMPDIR="/dev/shm/$dataset_id"
    fi
    if [ ! -d \$TMPDIR ]; then mkdir -p \$TMPDIR; fi

    STAR --genomeDir ${genomeIndex} \
        --readFilesIn $r1 $r2 \
        --readFilesCommand zcat \
        --sjdbGTFfile ${gtf} \
        --runThreadN ${task.cpus} \
        --twopassMode Basic \
        --outSAMunmapped Within \
        --outFileNamePrefix "${dataset_id}." \
        --chimSegmentMin 12 \
        --chimJunctionOverhangMin 12 \
        --alignSJDBoverhangMin 10 \
        --alignMatesGapMax 100000 \
        --chimSegmentReadGapMax 3 \
        --outFilterMismatchNmax 2 \
        --outSAMtype BAM Unsorted \
        --quantMode TranscriptomeSAM 
    samtools sort -@ ${task.cpus} -T \$TMPDIR -o ${dataset_id}.Aligned.sortedByCoord.out.bam -O BAM ${dataset_id}.Aligned.out.bam
    samtools index -@ ${task.cpus} ${dataset_id}.Aligned.sortedByCoord.out.bam
    """
}
