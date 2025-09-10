<h1 align="center">🧬CCRGeneticsBranch / Oncogenomics_NF_WF</h1>
<p align="center">
  A containerized Nextflow pipeline for scalable oncogenomics analyses
</p>

<p align="center">
  <a href="https://www.nextflow.io/">
    <img src="https://img.shields.io/badge/Nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg?labelColor=2b2d42&logo=nextflow&logoColor=white" alt="Nextflow">
  </a>
  <a href="https://www.docker.com/">
    <img src="https://img.shields.io/badge/Run%20with-Docker-2496ED.svg?labelColor=2b2d42&logo=docker&logoColor=white" alt="Docker">
  </a>
  <a href="https://sylabs.io/docs/">
    <img src="https://img.shields.io/badge/Run%20with-Singularity-1d355c.svg?labelColor=2b2d42&logo=singularity&logoColor=white" alt="Singularity">
  </a>
  <a href="https://github.com/CCRGeneticsBranch/Oncogenomics_NF_WF/releases/latest">
    <img src="https://img.shields.io/github/v/release/CCRGeneticsBranch/Oncogenomics_NF_WF?color=ff7b00&label=Latest%20Release&labelColor=2b2d42&logo=github" alt="Latest Release">
  </a>
  <a href="https://ccrgeneticsbranch.github.io/Oncogenomics_NF_WF/">
    <img src="https://img.shields.io/badge/Docs-Oncogenomics_NF_WF-008080.svg?labelColor=2b2d42&logo=readthedocs&logoColor=white" alt="Documentation">
  </a>
</p>

The Oncogenomics_NF_WF is a containerized Nextflow pipeline for processing exome and RNA-seq cancer data. It is built for scalable execution on HPC (Biowulf) and AWS. It integrates tools for variant calling, CNV detection, mutational signatures, TMB, HLA typing, neoantigen prediction, RNA quantification, fusion detection, and immune infiltration metrics. It currently supports hg19 and mm39 reference genomes. Efforts ongoing to expand support for hg38. For GRCm39 (mm39) reference, only mapping and RSEM quantification is supported.

## Requirements

- Nextflow (>= 23.10.0)
- Graphviz 2.40
- Docker or Singularity 3.10.5
- Access to Biowulf HPC or AWS Batch

## Usage

```
# Run on Biowulf
# Copy the launcher from the Khanlab space and run with a samplesheet
# cp /data/khanlab/projects/Nextflow_dev/dev/AWS_POC_MVP_NF/launch.sh .
./launch.sh /path/to/samplesheet.csv
```

## Input / Output

Input: Sample sheet with paths to fastq or bam/cram files. Link to [samplesheet guidelines](https://ccrgeneticsbranch.github.io/Oncogenomics_NF_WF/usage/).

Output: VCFs, CNV calls, mutational signatures, TMB metrics, HLA predictions, neoantigen candidates, expression matrices, fusion calls, immune infiltration metrics, QC reports.

Results are organized by patientID/casename with an option to visualize data at [clinOmics data portal](https://oncogenomics.ccr.cancer.gov/production/public/)

## Documentation

Full documentation is available here:
👉 [Oncogenomics_NF_WF Documentation](https://ccrgeneticsbranch.github.io/Oncogenomics_NF_WF/)

## Contributing

Issues and pull requests are welcome. Please follow coding and documentation guidelines.

Please send your comments/questions/suggestions to [Vineela Gangalapudi](https://github.com/vinegang) via [email](mailto:vineela.gangalapudi@nih.gov).
