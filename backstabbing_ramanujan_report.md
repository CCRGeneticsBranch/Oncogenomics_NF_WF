## cutadapt (Test2_R_T)

script:

    
	cutadapt --pair-filter=any \
	--nextseq-trim=2 \
	--trim-n \
	-n 5 -O 5 \
	-q 10,10 -m 30:30 \
	-b file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa \
	-B file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa \
	-j 4 \
        -o Test2_R_T_R1.trim.fastq.gz \
        -p Test2_R_T_R2.trim.fastq.gz \
	Test2_R_T_R1.fastq.gz Test2_R_T_R2.fastq.gz      
        

exit status: 0
task status: COMPLETED
task hash: ef/130a8d
task name: cutadapt (Test2_R_T)

## cutadapt (Test3_R_T)

script:

    
	cutadapt --pair-filter=any \
	--nextseq-trim=2 \
	--trim-n \
	-n 5 -O 5 \
	-q 10,10 -m 30:30 \
	-b file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa \
	-B file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa \
	-j 4 \
        -o Test3_R_T_R1.trim.fastq.gz \
        -p Test3_R_T_R2.trim.fastq.gz \
	Test3_R_T_R1.fastq.gz Test3_R_T_R2.fastq.gz      
        

exit status: 0
task status: COMPLETED
task hash: e8/8bf60f
task name: cutadapt (Test3_R_T)

## cutadapt (Test1_R_T)

script:

    
	cutadapt --pair-filter=any \
	--nextseq-trim=2 \
	--trim-n \
	-n 5 -O 5 \
	-q 10,10 -m 30:30 \
	-b file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa \
	-B file:/opt2/TruSeq_and_nextera_adapters.consolidated.fa \
	-j 4 \
        -o Test1_R_T_R1.trim.fastq.gz \
        -p Test1_R_T_R2.trim.fastq.gz \
	Test1_R_T_R1.fastq.gz Test1_R_T_R2.fastq.gz      
        

exit status: 0
task status: COMPLETED
task hash: a6/09e472
task name: cutadapt (Test1_R_T)

