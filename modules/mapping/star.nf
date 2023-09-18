process Star {
    tag "$meta.lib"
    scratch true
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}", mode: "copy"

    input:
    tuple val(meta), path(trim),path(star_genomeIndex),path(gtf)
        
    output:
    tuple val(meta), path("${meta.lib}.Aligned.toTranscriptome.out.bam") , emit: transcriptome_bam
    tuple val(meta), path("${meta.lib}.Aligned.sortedByCoord.out.bam"), emit: genome_bam
    tuple val(meta), path("${meta.lib}.Aligned.sortedByCoord.out.bam.bai"), emit: genome_bai
    tuple val(meta), path("${meta.lib}.Chimeric.out.junction"), emit: chimeric_junction

    stub:
    """
    touch "${meta.lib}.Aligned.toTranscriptome.out.bam"
    touch "${meta.lib}.Aligned.sortedByCoord.out.bam"
    touch "${meta.lib}.Aligned.sortedByCoord.out.bam.bai"
    touch "${meta.lib}.Chimeric.out.junction"
    """
    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    set -exo pipefail
    echo $PWD
    echo $SLURM_JOB_ID
    
        
        # run STAR alignment
        STAR --genomeDir ${star_genomeIndex} \
            --readFilesIn ${trim[0]} ${trim[1]} \
            --readFilesCommand zcat \
            --sjdbGTFfile ${gtf} \
            --runThreadN ${task.cpus} \
            --twopassMode Basic \
            --outSAMunmapped Within \
            --outFileNamePrefix ${prefix}. \
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
        samtools sort -@ !{task.cpus}  -o ${prefix}.Aligned.sortedByCoord.out.bam -O BAM ${prefix}.Aligned.out.bam
    

    # index files
    samtools index -@ ${task.cpus} ${prefix}.Aligned.sortedByCoord.out.bam
    """
}
