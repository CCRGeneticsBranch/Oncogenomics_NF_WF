process GATK_RNASeq_Trim {
        tag { dataset_id }

        input:
        tuple val(dataset_id),
        path(bam),
        path(index),
        path(markdup),
        path(genome),
        path(genome_fai),
        path(genome_dict)

        output:
        tuple val("${dataset_id}"),
        path("trim_${dataset_id}.star.trim.bam"),
        path("trim_${dataset_id}.star.trim.bai")


        script:
        """

        java -jar /opt2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T SplitNCigarReads -R $genome -I $bam -o trim_${dataset_id}.star.trim.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

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

	java -jar /opt2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T RealignerTargetCreator -nt 10 -R $genome -known $phase1_1000g -known $Mills_and_1000g -I $bam -o trim_${dataset_id}.star.realignment.intervals

        java -jar /opt2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T IndelRealigner -R $genome -known $phase1_1000g -known $Mills_and_1000g -I $bam --targetIntervals trim_${dataset_id}.star.realignment.intervals -o trim_${dataset_id}.star.Ir.bam        

        """
}



process GATK_RNASeq_BR_PR {
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
        path("trim_${dataset_id}.star.final.bam"),
        path("trim_${dataset_id}.star.final.bam.bai")


        script:
        """

        java -jar /opt2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T BaseRecalibrator -R $genome -knownSites $phase1_1000g -knownSites $Mills_and_1000g -I $bam -o trim_${dataset_id}.star.recalibration.matrix.txt 

        java -jar /opt2/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T PrintReads -R $genome -I $bam -o trim_${dataset_id}.star.final.bam -BQSR trim_${dataset_id}.star.recalibration.matrix.txt

        mv trim_${dataset_id}.star.final.bai trim_${dataset_id}.star.final.bam.bai

        """
}
