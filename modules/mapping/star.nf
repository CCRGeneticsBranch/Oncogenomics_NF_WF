process Star {
    tag "$meta.lib"
    scratch true
    publishDir "${params.resultsdir}/${meta.id}/${meta.casename}/${meta.lib}", mode: "${params.publishDirMode}" ,pattern: "${meta.lib}*bed.gz*"

    input:
    tuple val(meta), path(trim),path(star_genomeIndex),path(gtf)

    output:
    tuple val(meta), path("${meta.lib}.Aligned.toTranscriptome.out.bam") , emit: transcriptome_bam
    tuple val(meta), path("${meta.lib}.Aligned.sortedByCoord.out.bam"), emit: genome_bam
    tuple val(meta), path("${meta.lib}.Aligned.sortedByCoord.out.bam.bai"), emit: genome_bai
    tuple val(meta), path("${meta.lib}.Chimeric.out.junction"), emit: chimeric_junction
    tuple val(meta), path("${meta.lib}.SJ.out.bed.gz"), emit: sj_bed
    tuple val(meta), path("${meta.lib}.SJ.out.bed.gz.tbi"), emit: sj_bed_tbi
    path "versions.yml"             , emit: versions

    stub:
    """
    touch "${meta.lib}.Aligned.toTranscriptome.out.bam"
    touch "${meta.lib}.Aligned.sortedByCoord.out.bam"
    touch "${meta.lib}.Aligned.sortedByCoord.out.bam.bai"
    touch "${meta.lib}.Chimeric.out.junction"
    touch "${meta.lib}.SJ.out.bed.gz"
    touch "${meta.lib}.SJ.out.bed.gz.tbi"
    """
    script:
    def prefix = task.ext.prefix ?: "${meta.lib}"
    """
    TMP=tmp/
    mkdir -p \$TMP
    trap 'rm -rf "\$TMP"' EXIT


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
        samtools sort -@ ${task.cpus}  -T \$TMP -o ${prefix}.Aligned.sortedByCoord.out.bam -O BAM ${prefix}.Aligned.out.bam


        awk -F"\t" 'BEGIN{OFS="\t"}{strand=(\$4==1)?"+":"-";annotated=(\$6==1)?"true":"false";if(\$5==0) motif="non-canonical"; if(\$5==1)motif="GT/AG";if(\$5==2)motif="CT/AC";if(\$5==3)motif="GC/AC";if(\$5==4)motif="CT/GC";if(\$5==5)motif="AT/AC";if(\$5==6)motif="GT/AT";print \$1,\$2,\$3,"motif="motif";uniquely_mapped="\$7";multi_mapped="\$8";maximum_spliced_alignment_overhang="\$9";annotated_junction="annotated,\$7,strand}' \
        ${prefix}.SJ.out.tab | bedtools sort -i - | bgzip > ${prefix}.SJ.out.bed.gz
        tabix -0 -p bed ${prefix}.SJ.out.bed.gz

    # index files
    samtools index -@ ${task.cpus} ${prefix}.Aligned.sortedByCoord.out.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        STAR: \$(STAR --version)
    END_VERSIONS
    """
}
