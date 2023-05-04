process GATK_RNASeq_Trim {
        tag { dataset_id }

        input:
        tuple val(dataset_id),
        val(library),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict)

        output:
        tuple val("${dataset_id}"),
        val("$library"),
        path("${library}.star.trim.bam"),
        path("${library}.star.trim.bai")


        script:
        """

        java -jar -Xmx40g \$GATK_JAR -T SplitNCigarReads -R $genome -I $bam -o ${library}.star.trim.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

        """
}


process GATK_RNASeq_RTC_IR {
        tag { dataset_id }

        input:
	tuple val(dataset_id),
        val(library),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(phase1_1000g),
        path(Mills_and_1000g)

        output:
        tuple val("${dataset_id}"),
        val("$library"),
        path("${library}.star.Ir.bam"),
        path("${library}.star.Ir.bai")


        script:
        """

	java -jar -Xmx40g \$GATK_JAR -T RealignerTargetCreator -nt 10 -R $genome -known $phase1_1000g -known $Mills_and_1000g -I $bam -o ${library}.star.realignment.intervals

        java -jar -Xmx40g \$GATK_JAR -T IndelRealigner -R $genome -known $phase1_1000g -known $Mills_and_1000g -I $bam --targetIntervals ${library}.star.realignment.intervals -o ${library}.star.Ir.bam        

        """
}



process GATK_RNASeq_BR_PR {
        tag { dataset_id }

        publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/${library}", mode: "${params.publishDirMode}"

        input:
	tuple val(dataset_id),
        val(library),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(phase1_1000g),
        path(Mills_and_1000g)

        output:
        tuple val("${dataset_id}"),
        val("$library"),
        path("${library}.star.final.bam"),
        path("${library}.star.final.bam.bai")


        script:
        """

        java -jar -Xmx40g \$GATK_JAR -T BaseRecalibrator -R $genome -knownSites $phase1_1000g -knownSites $Mills_and_1000g -I $bam -o ${library}.star.recalibration.matrix.txt 

        java -jar -Xmx40g \$GATK_JAR -T PrintReads -R $genome -I $bam -o ${library}.star.final.bam -BQSR ${library}.star.recalibration.matrix.txt

        mv ${library}.star.final.bai ${library}.star.final.bam.bai

        """
}


process RNAseq_HaplotypeCaller {

     tag { dataset_id }

     publishDir "${params.resultsdir}/${dataset_id}/${params.casename}/${library}/calls", mode: "${params.publishDirMode}"

     input:
     tuple val(dataset_id),
        val(library),
        path(bam),
        path(index),
        path(genome),
        path(genome_fai),
        path(genome_dict),
        path(dbsnp)

     output:
     tuple val("${dataset_id}"),
        val("$library"),
        path("${library}.HC_RNASeq.raw.vcf")

     stub:
     """
     touch "${library}.HC_RNASeq.raw.vcf"
     """

     shell:
     '''
     set -exo pipefail
     java -jar -Xmx40g \$GATK_JAR -T HaplotypeCaller -R !{genome} -I !{bam} -o !{library}.vcf --dbsnp !{dbsnp} -dontUseSoftClippedBases -stand_call_conf 30 -nct !{task.cpus}
     java -jar -Xmx40g \$GATK_JAR -T VariantFiltration -R !{genome} -V !{library}.vcf -window 35 -cluster 3 --filterExpression "FS > 30.0 || QD < 2.0" -filterName "RNASeqFilters_FS_QD" -o !{library}.HC_RNASeq.raw.vcf
     '''
}


