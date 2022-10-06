params.reads = "s3://ccr-genomics-testdata/testdata/Test*_R_T_R{1,2}.fastq.gz"
reads_ch = Channel.fromFilePairs(params.reads)

process cutadapt {
	container 'nciccbr/ncigb_cutadapt_v1.18:latest'
	publishDir 's3://agc-424336837382-us-east-1/nfmvpout', mode: 'copy'
	
	input:
	tuple val(sample_id), path(reads) from reads_ch

	output:
	path "trim*" into trim_ch1, trim_ch2

	
	"""
	cutadapt  -o trim_${sample_id}_R1.fastq -p trim_${sample_id}_R2.fastq ${reads[0]} ${reads[1]}
	"""

}
