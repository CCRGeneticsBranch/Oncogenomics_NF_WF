## Results

Results from the pipeline are organized in a structured directory format, making it easy to locate and analyze outputs.

### Results Directory Structure

The results root is defined by `${OUTDIR}`, which can be either the default location or a custom path specified when launching the workflow. The directory structure is as follows:

    ├── Actionable
    ├── annotation
    ├── log
    ├── Test
    ├── Test_N2D_E
    ├── Test_T1D_E
    ├── Test_T1R_T
    ├── qc
    └── successful.txt

**Description of Contents**


**Actionable/**  
 
??? "Contains Actionable fusion calls and mutational signatures from COSMIC"

    - `*.unionSomaticVarsFull.txt` – Union of somatic variants from multiple variant callers.  
    - `*.mutationalSignature.pdf` – Mutational signature plots derived from COSMIC signatures.  
    - `*.actionable.fusion.txt` – Actionable fusion calls from consensus fusion detection.  
    - `fusion.actionable.txt` – Summary of actionable fusion results across samples.


 **annotation/**  
??? "Contains variant annotations from multiple databases and tools:"

    - `AnnotationInput` – Base variant list before annotation.  
    - `AnnotationInput.anno` – Primary annotation file with basic functional effects.  
    - `AnnotationInput.sift` – SIFT scores predicting functional impact.  
    - `AnnotationInput.candl` – CAN-DL annotation results.  
    - `AnnotationInput.civic` – CIViC (Clinical Interpretation of Variants in Cancer) annotations.  
    - `AnnotationInput.clinvar` – ClinVar clinical significance annotations.  
    - `AnnotationInput.docm` – Database of Curated Mutations results.  
    - `AnnotationInput_final` – Consolidated annotation file after merging data sources.  
    - `AnnotationInput.hgmd` – Human Gene Mutation Database annotations.  
    - `AnnotationInput.match` – Matched database entries.  
    - `AnnotationInput.mcg` – MCG (Molecular Cancer Genetics) annotations.  
    - `AnnotationInput.tcc` – TCC-specific annotations.  
    - `AnnotationInput.cadd` – CADD (Combined Annotation Dependent Depletion) scores.  
    - `AnnotationInput.clinseq` – ClinSeq annotations.  
    - `AnnotationInput.cosmic` – COSMIC (Catalogue of Somatic Mutations in Cancer) annotations.  
    - `AnnotationInput.gene` – Gene-level annotations.  
    - `AnnotationInput.pcg` – PCG (Precision Cancer Genomics) annotations.  
    - `versions.yml` – Software and database version information.  
    - `*.Annotations.coding.rare.txt` – Coding variants that are rare in the population with MAF < 0.05. 
    - `*.Annotations.final.txt` – Final comprehensive variant annotation file.

**log/**

??? "Nextflow run logs and reports"

    - `trace.YYYY-MM-DD_HH-MM-SS.txt` – Nextflow trace file showing per-process execution stats (CPU, memory, time, etc.).
    - `report_YYYY-MM-DD_HH-MM-SS.html` – Interactive Nextflow HTML report summarizing execution, resources, and task breakdowns.
    - `timeline_YYYY-MM-DD_HH-MM-SS.html` – Visual timeline showing when each process ran and for how long.
    - `dag.YYYY-MM-DD_HH-MM-SS.png` – Directed acyclic graph of the workflow steps and their dependencies.
    - `bco_YYYY-MM-DD_HH-MM-SS.json` – BioCompute Object (BCO) file capturing workflow provenance and metadata.
    - `manifest_YYYY-MM-DD_HH-MM-SS.json` – Manifest describing run inputs, parameters, and environment details.

**Test_N2D_E/**


??? "Matched normal library outputs"

    - **calls/** – Variant call files generated for the matched normal sample.
    - **HLA/** – HLA typing results.
    - `Test_N2D_E.final.bam` – Final aligned GATK BAM file for the matched normal sample.
    - `Test_N2D_E.final.bam.bai` – Index file for the final BAM.
    - `Test_N2D_E.final.squeeze.bam` – BAM for IGV.
    - `Test_N2D_E.final.squeeze.bam.bai` – Index file for the squeezed BAM.
    - **qc/** – Quality control metrics and reports.
    - **TCellExTRECT/** – T cell infiltration/extrect analysis results.
    - **verifyBamID/** – Results from VerifyBamID for contamination/identity checks.
  
**Test_T1D_E/**


??? "Tumor DNA library outputs"

    - **calls/** – Variant call files for the tumor DNA sample.
    - **cnvkit/** – CNVkit copy-number analysis results.
    - **HLA/** – HLA typing results.
    - `Test_T1D_E.final.bam` – Final aligned GATK BAM file for the tumor DNA sample.
    - `Test_T1D_E.final.bam.bai` – Index file for the final BAM.
    - `Test_T1D_E.final.squeeze.bam` – BAM for IGV viewing.
    - `Test_T1D_E.final.squeeze.bam.bai` – Index file for the squeezed BAM.
    - **NeoAntigen/** – Predicted neoantigen results.
    - **qc/** – Quality control metrics and reports.
    - **sequenza/** – Sequenza purity/ploidy estimation results.
    - **TCellExTRECT/** – T cell infiltration/extrect analysis results.
    - **verifyBamID/** – VerifyBamID contamination/identity check results.


**Test_T1R_T/**

??? "Tumor RNA library outputs"
    - **calls/** – Variant calls from RNA-seq data.
    - **fusion/** – Fusion detection results.
    - **HLA/** – HLA typing results from RNA data.
    - `Test_T1R_T.final.bam` / `.bai` – Final aligned RNA BAM + index.
    - `Test_T1R_T.final.squeeze.bam` / `.bai` – Squeezed BAM for IGV + index.
    - **qc/** – QC metrics/reports for RNA-seq.
    - **RSEM_ENS/** – RSEM expression quantification.
  
**qc/**


??? "Patient-level QC results (click to expand)"
    - `Test.transcriptCoverage.png` – Transcript coverage plot for RNA-seq libraries.
    - `Test.RnaSeqQC.txt` – RNA-seq QC metrics table.
    - `Test.circos.png` – Circos plot summarizing structural variations and genome-wide metrics across libraries within the case.
    - `Test.hotspot_coverage.png` – Coverage at hotspot variant locations across libraries within the case.
    - `Test.consolidated_QC.txt` – Consolidated QC summary combining results from all libraries for the patient.
    - `multiqc_report.html` – Aggregated QC metrics from multiple tools in HTML format.
    - `Test.config.*.txt` – Tool versions used for each run.
    - `genotyping.html` – Genotyping QC results in HTML format to send out pipeline completion emails.
    - `Test.genotyping.txt` – Genotyping QC results in text format.

**successful.txt** - This file acts as a trigger to upload results to the [ClinOmics Data Portal](https://oncogenomics.ccr.cancer.gov/production/public/).