process GATK_RNASeq_Trim {
        tag { dataset_id }

        input:
        tuple val(dataset_id),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict)

        output:
        tuple val("${dataset_id}"),
        path("trim_${dataset_id}.star.trim.bam"),
        path("trim_${dataset_id}.star.trim.bai")


        script:
        """

        java -jar \$GATK_JAR -T SplitNCigarReads -R $genome -I $bam -o trim_${dataset_id}.star.trim.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

        """
}


process GATK_RNASeq_RTC_IR {
        tag { dataset_id }

        input:
	tuple val(dataset_id),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(phase1_1000g),
        path(Mills_and_1000g)

        output:
        tuple val("${dataset_id}"),
        path("trim_${dataset_id}.star.Ir.bam"),
        path("trim_${dataset_id}.star.Ir.bai")


        script:
        """

	java -jar \$GATK_JAR -T RealignerTargetCreator -nt 10 -R $genome -known $phase1_1000g -known $Mills_and_1000g -I $bam -o trim_${dataset_id}.star.realignment.intervals

        java -jar \$GATK_JAR -T IndelRealigner -R $genome -known $phase1_1000g -known $Mills_and_1000g -I $bam --targetIntervals trim_${dataset_id}.star.realignment.intervals -o trim_${dataset_id}.star.Ir.bam        

        """
}



process GATK_RNASeq_BR_PR {
        tag { dataset_id }

        publishDir "${params.resultsdir}/${dataset_id}/GATK_RNAseq", mode: "${params.publishDirMode}"

        input:
	tuple val(dataset_id),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(phase1_1000g),
        path(Mills_and_1000g)

        output:
        tuple val("${dataset_id}"),
        path("trim_${dataset_id}.star.final.bam"),
        path("trim_${dataset_id}.star.final.bam.bai")


        script:
        """

        java -jar \$GATK_JAR -T BaseRecalibrator -R $genome -knownSites $phase1_1000g -knownSites $Mills_and_1000g -I $bam -o trim_${dataset_id}.star.recalibration.matrix.txt 

        java -jar \$GATK_JAR -T PrintReads -R $genome -I $bam -o trim_${dataset_id}.star.final.bam -BQSR trim_${dataset_id}.star.recalibration.matrix.txt

        mv trim_${dataset_id}.star.final.bai trim_${dataset_id}.star.final.bam.bai

        """
}


process RNAseq_HaplotypeCaller {

     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/GATK_RNAseq", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(dbsnp)

     output:
     tuple val("${dataset_id}"),
        path("${dataset_id}.HC_RNASeq.raw.vcf")

     stub:
     """
     touch "${dataset_id}.HC_RNASeq.raw.vcf"
     """

     shell:
     '''
     set -exo pipefail
     java -jar \$GATK_JAR -T HaplotypeCaller -R !{genome} -I !{bam} -o !{dataset_id}.vcf --dbsnp !{dbsnp} -dontUseSoftClippedBases -stand_call_conf 30 -nct !{task.cpus}
     java -jar \$GATK_JAR -T VariantFiltration -R !{genome} -V !{dataset_id}.vcf -window 35 -cluster 3 --filterExpression "FS > 30.0 || QD < 2.0" -filterName "RNASeqFilters_FS_QD" -o !{dataset_id}.HC_RNASeq.raw.vcf
     '''
}


