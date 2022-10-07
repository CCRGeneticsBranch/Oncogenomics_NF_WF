// Using DSL-2
nextflow.enable.dsl=2

params.reads = "s3://ccr-genomics-testdata/testdata/Test*_R_T_R{1,2}.fastq.gz"
//reads_ch = Channel.fromFilePairs(params.reads)

//Print out log
log.info """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         reads        : ${params.reads}
         """
         .stripIndent()

// set up processes
process cutadapt {
        tag { dataset_id }
	publishDir 's3://agc-424336837382-us-east-1/nfmvpout', mode: 'copy'

        input:
        tuple val(dataset_id), 
        path(forward), 
        path(reverse)

        output:
        tuple val("${dataset_id}"),
        path("trim_${dataset_id}_R1.fastq"),
        path("trim_${dataset_id}_R2.fastq")

        container 'nciccbr/ncigb_cutadapt_v1.18:latest'

        script:
        """
	cutadapt  -o trim_${dataset_id}_R1.fastq -p trim_${dataset_id}_R2.fastq $forward $reverse
        """

}

process fastqc {
        tag { dataset_id }
        publishDir 's3://agc-424336837382-us-east-1/nfmvpout', mode: 'copy'

//	cache false
        input:
        tuple val(dataset_id),
        path(forward),
        path(reverse)

        output:
        tuple val("${dataset_id}"),
        path("fastqc_trim_${dataset_id}")

        container 'nciccbr/ccbr_fastqc_0.11.9:v1.1'

        script:
        """	
        mkdir fastqc_trim_${dataset_id}
        fastqc -o fastqc_trim_${dataset_id} -q $forward $reverse
        """
}

workflow{
    read_pairs = Channel
        .fromFilePairs(params.reads, flat: true)
        .ifEmpty { exit 1, "Read pairs could not be found: ${params.reads}" }

    cutadapt(read_pairs)
    fastqc(cutadapt.out)
}

