# Early-life Iron Supplementation Impact on the Gut Microbiota  
*A pipeline for microbiota analysis*

<p align="center">
  <img src="https://github.com/bioth/gut-microbiota-iron/blob/scripts/photos/pipeline.png?raw=true"/>
</p>

This is a walkthrough of the analysis pipeline built for analysis of the NovaSeq sequencing of the 16S V5–V6 region. It covers preprocessing, ASV construction, taxonomic annotation, and downstream visualizations and analyses.

---

## Requirements

Main packages and softwares required to run the pipeline include:

- **R** (version ≥ 4.3.0)  
- **R tools**  
- **BiocManager** (v3.18)  
- **DADA2** (v3.18)  
- **ggplot2** (v3.4.4)  
- **phyloseq** (v1.52)  
- **DESeq2** (v1.48.1)  
- **Cutadapt** (v4.8)  

---

## 1. Primer Trimming

Primer sequences were removed with:

[cutAdapt.sh](1-Primer%20trimming/cutAdapt.sh)

This script leverages Cutadapt in paired-end mode to trim primers efficiently.

---

## 2. DADA2 – ASV Construction & Taxonomic Annotation

ASV inference—including filtering, denoising, and chimera removal—was executed with:

[1-Building ASVs dataset.R](2-DADA/1-Building%20ASVs%20dataset.R)

A custom error model was implemented to account for binned quality scores in FASTQ files, enabling more accurate error modeling ([DADA2 issue #1307](https://github.com/benjjneb/dada2/issues/1307)).

Execution occurred on Alliance Canada’s HPC infrastructure to leverage its computational power.

### Quality Profiles
<p align="center">
  <img src="https://github.com/bioth/gut-microbiota-iron/blob/scripts/photos/quality_profileR1.png?raw=true" width="45%">
  <img src="https://github.com/bioth/gut-microbiota-iron/blob/scripts/photos/quality_profileR2.png?raw=true" width="45%">
</p>

*Figure: Quality profiles of forward (R1) and reverse (R2) reads. The plots show per-base quality scores across sequencing cycles, used to select truncation and filtering parameters in DADA2.*

### Error Model
<p align="center">
  <img src="https://github.com/bioth/gut-microbiota-iron/blob/scripts/photos/dada_plotR1_m1.png?raw=true" width="45%">
  <img src="https://github.com/bioth/gut-microbiota-iron/blob/scripts/photos/dada_plotR2_m1.png?raw=true" width="45%">
</p>

*Figure: DADA2 error models for forward (R1) and reverse (R2) reads. Black points show observed error rates for each nucleotide transition, while the red line represents the fitted error model used during ASV inference.*

Taxonomic annotation was performed via:

[2-Taxonomic annotation.R](2-DADA/2-Taxonomic%20annotation.R)

The script uses the **SILVA v138.1** database. See [Silva 138.1 prokaryotic SSU taxonomic training data formatted for DADA2](https://zenodo.org/records/4587955).

---

## 3. Downstream Analysis

Alpha diversity, beta diversity, and differential abundance analyses were run with:

[microbiota_analysis.R](3-Analysis/microbiota_analysis.R)

Custom functions for visualizations, stats, and analysis can be found in the [`utils`](utils) folder.

---

## Extensions of a previously existing package for microbiota visualization

###  StackbarExtended

We leveraged the [**StackbarExtended**](https://github.com/ThibaultCuisiniere/StackbarExtended) R package—a peer-reviewed tool for visualizing taxa relative abundance with integrated phylogenetic information and differential abundance statistics.

- We extended its main functionality to accommodate our experimental design (~ treatment * diet).
- The customization script is available here: [plot_microbiota_extension.R](utils/plot_microbiota_extension.R).

<p align="center">
  <img src="https://github.com/bioth/gut-microbiota-iron/blob/scripts/photos/stackbar.png?raw=true" width="70%">
</p>

---

###  Chronobiome (Custom tool)

Inspired by StackbarExtended’s logic to assign colors to taxa, we developed a custom framework—**Chronobiome**—for visualizing longitudinal microbiota data. Features include:

- Group-averaged relative abundance over time  
- Flexible taxonomic-depth views (either comprehensive or targeted)  

<p align="center">
  <img src="https://github.com/bioth/gut-microbiota-iron/blob/scripts/photos/chronobiome.png?raw=true" width="70%">
</p>

<p align="center">
  <img src="https://github.com/bioth/gut-microbiota-iron/blob/scripts/photos/e.png?raw=true" width="70%">
</p>
