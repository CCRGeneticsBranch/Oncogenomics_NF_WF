## References

This pipeline supports the human genome **hg19 (GRCh37)** and the mouse genome **mm39 (GRCm39)**, using reference builds and annotations obtained from **GENCODE**.

### hg19 (GRCh37) resources

| Component | File | Source |
|---|---|---|
| Reference genome (FASTA) | `GRCh37.primary_assembly.genome.fa.gz` | ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh37_mapping/GRCh37.primary_assembly.genome.fa.gz |
| Gene annotation (GTF) | `gencode.v37lift37.annotation.gtf.gz` | ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh37_mapping/gencode.v37lift37.annotation.gtf.gz |

> ERCC spike-in sequences were added to this reference (FASTA concatenated and indices rebuilt accordingly).

### mm39 (GRCm39) resources

| Component | File | Source |
|---|---|---|
| Reference genome (FASTA) | `GRCm39.primary_assembly.genome.fa.gz` | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/GRCm39.primary_assembly.genome.fa.gz |
| Gene annotation (GTF) | `gencode.vM33.primary_assembly.annotation.gtf.gz` | https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.primary_assembly.annotation.gtf.gz |




### Tools Used

| Tool | Version |
|---|---|
| Arriba | 2.3.0 |
| bwa | 0.7.17-r1188 |
| samtools | 1.15.1 |
| cnvkit | 0.9.6 |
| bedtools | v2.27.1 |
| Cutadapt | 1.18 |
| GATK | 3.8-1-0-gf15c1c3ef |
| Fastqc | v0.11.9 |
| Fusioncatcher | 1.33 |
| GATK | 3.8-1-0-gf15c1c3ef |
| bcftools | 1.10.2 |
| HLA-HD | 1.7.0 |
| Picard | 2.27.4 |
| Kraken | 1.1.1 |
| manta | 1.6.0 |
| multiqc | 1.9 |
| mutect | 3.1-0-g72492bb |
| Nextflow | 23.10.0 |
| Optitype | v1.3.1 |
| GATK | 3.8-1-0-gf15c1c3ef |
| samtools | 1.10 |
| RSEM | v1.3.1 |
| sequenza-utils | 3.0.0 |
| Snpsift | 4.3t |
| python | 3.9.5 |
| yaml | 5.4.1 |
| STAR | 2.7.10a |
| STAR-Fusion | 1.11.1 |
| strelka | 2.9.10 |
| vcftools | 0.1.16 |
| VEP | 102.0 |
| verifyBamID | 1.1.3 |
