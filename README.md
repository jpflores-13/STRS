# Hyperosmotic stress induces chromatin loop reorganization through CTCF retention

This repository contains code, analysis scripts, and figure generation for our study investigating how hyperosmotic stress affects 3D genome organization and chromatin loop dynamics.

> **Flores JP**, [Additional Authors TBD]  
> *Hyperosmotic stress induces chromatin loop reorganization through CTCF retention*  
> Preprint: **TBD DOI** · Journal: **TBD**

## Overview

We integrated Hi-C, CUT&Tag, and RNA-seq data to characterize how sorbitol treatment (hyperosmotic stress) affects chromatin loop formation and gene expression in HEK293 cells. Our analysis reveals large-scale rewiring of chromatin interactions with distinct patterns of loop gain and loss, driven by differential CTCF and cohesin (RAD21) retention at loop anchors.

**Key findings:**
- Identification of gained and lost chromatin loops under hyperosmotic stress
- CTCF retention at loop anchors correlates with loop stability
- YAP1 and H3K27ac enrichment at gained loop anchors
- Timecourse analysis reveals rapid loop reorganization kinetics
- Cross-cell line validation in HEK293 WT, HCT116, and T47D cells

---

## Repository Structure

```
STRS/
├── scripts/
│   ├── analysis/         # Differential analysis, enrichment, statistics
│   ├── figures/          # Main and supplementary figure generation (R)
│   ├── processing/       # Loop extraction, DESeq2, APA calculations
│   └── utils/            # Helper functions and themes
├── data/                 # Raw and processed data (see GEO)
│   ├── raw/             # FASTQ, loop calls, Hi-C maps (not in repo)
│   └── processed/       # Differential loops, APA matrices (not in repo)
├── figures/             # Final manuscript figures (PDF)
├── tables/              # Gene lists and summary tables
├── plots/               # Exploratory plots (not in repo)
├── Makefile             # Figure generation automation
├── renv.lock            # R package versions for reproducibility
└── README.md            # This file
```

> **Note:** Large data files are not versioned in this repository. Raw and processed data are available via **GEO: GSE_TBD**. See `data/*/README.md` for details.

---

## Quick Start

### Prerequisites

- R ≥ 4.3.0
- Required R packages (managed via `renv`)
- GNU Make (optional, for automated figure generation)

### Setup

```bash
# Clone the repository
git clone https://github.com/YOUR_USERNAME/STRS.git
cd STRS

# Restore R package environment
R -e "install.packages('renv'); renv::restore()"

# Download processed data from GEO (once available)
# See data/README.md for instructions

# Generate all manuscript figures
make
```

### Generate Figures

```bash
make              # Generate all figures
make main         # Main figures (1-5) only
make supplementary # Supplementary figures only
make fig1         # Generate Figure 1 only
make clean        # Remove generated figures
make help         # See all available commands
```

---

## Data Availability

### Raw Data
- **Hi-C:** FASTQ files on SRA: **SRP_TBD**
- **CUT&Tag:** FASTQ files on SRA: **SRP_TBD**
- **RNA-seq:** FASTQ files on SRA: **SRP_TBD**

### Processed Data
All processed data available via **GEO: GSE_TBD**:
- Hi-C contact maps (`.hic` format)
- Loop calls (BEDPE format)
- Differential loop analysis results (`.rds`)
- Normalized APA matrices (`.rds`)
- CUT&Tag peak calls (narrowPeak format)
- CUT&Tag signal tracks (bigWig format)
- RNA-seq count matrices (`.rds`)
- DESeq2 results (`.rds`)

### Sample Metadata
- Sample sheets: `data/*/README.md`
- Library information: See GEO submission

---

## Computational Methods

### Hi-C Processing
- Loop calling: SIP (Significant Interaction Peak caller) via custom pipeline
- Differential loops: DESeq2 with replicate-aware design
- APA (Aggregate Peak Analysis): `mariner` R package
- Resolution: 10kb bins
- Normalization: Variance-stabilizing transformation for PCA; none for DESeq2

### CUT&Tag Analysis
- Peak calling: SEACR (Sparse Enrichment Analysis for CUT&RUN)
- Differential analysis: CPM binning strategy to avoid power bias
- Proteins analyzed: CTCF, RAD21, YAP1, H3K27ac
- Motif analysis: HOMER

### RNA-seq Analysis
- Alignment: STAR
- Quantification: featureCounts
- Differential expression: DESeq2 with likelihood ratio test (LRT)
- Clustering: K-means on variance-stabilized counts

### Software Versions
- R 4.3.1
- DESeq2 1.42.0
- mariner 1.2.0
- plotgardener 1.8.0
- See `renv.lock` for complete package versions

---

## Reproducibility

### R Environment Management
This project uses `renv` for reproducible package management:

```r
# Install exact package versions used in manuscript
renv::restore()

# Check package status
renv::status()

# Update lockfile after adding packages
renv::snapshot()
```

### Figure Generation
All figures can be regenerated from processed data:

```bash
# Individual figures
Rscript scripts/figures/Figure1.R
Rscript scripts/figures/Figure2.R
# ... etc

# Or use Make
make fig1 fig2 fig3
```

### Processing Pipeline
To reproduce processing from raw data (requires significant compute resources):

```r
# 1. Extract loop counts
source("scripts/processing/extractLoops.R")

# 2. Call differential loops
source("scripts/processing/callDiffLoops.R")

# 3. Calculate APA matrices
source("scripts/processing/calcAPA_gained.R")
source("scripts/processing/calcAPA_lost.R")
source("scripts/processing/calcAPA_timecourse_gained.R")
source("scripts/processing/calcAPA_timecourse_lost.R")
```

**Note:** Processing scripts require access to raw data and substantial computational resources (32+ GB RAM, several hours runtime). We recommend using processed data from GEO for figure reproduction.

---

## Manuscript Figures

### Main Figures
- **Figure 1:** Differential loop analysis and APA validation across cell lines
- **Figure 2:** Hi-C timecourse analysis of loop dynamics
- **Figure 3:** CUT&Tag enrichment at loop anchors (CTCF, RAD21, YAP1, H3K27ac)
- **Figure 4:** CTCF motif analysis and retention patterns
- **Figure 5:** RNA-seq timecourse and gene expression at loop anchors

### Supplementary Figures
- **Figure S1:** Loop overlap analysis and size distributions
- **Figure S2:** Additional Hi-C validation and controls
- **Figure S3:** Extended CUT&Tag and RNA-seq analyses

---

## Citation

If you use code or data from this repository, please cite:

```
Flores JP, et al. (2024)
Hyperosmotic stress induces chromatin loop reorganization through CTCF retention
[Journal TBD] [DOI TBD]
```

Additionally, please cite the software packages used:
- **mariner:** Kramer et al. (2022) https://doi.org/10.1093/bioinformatics/btac062
- **plotgardener:** Kramer et al. (2022) https://doi.org/10.1093/bioinformatics/btab761
- **DESeq2:** Love et al. (2014) https://doi.org/10.1186/s13059-014-0550-8

For a specific version of this repository, cite the Zenodo DOI: **TBD**

---

## Contributing

This is a research project repository associated with a manuscript. We welcome:
- **Bug reports:** Open an issue if you find errors in code or documentation
- **Questions:** Open an issue for clarification on methods or results
- **Suggestions:** Propose improvements via pull request

For major changes, please open an issue first to discuss proposed modifications.

---

## License

- **Code:** MIT License (or as specified in individual files)
- **Figures and data outputs:** CC BY-NC 4.0
- **Raw data:** See GEO submission for data use policies

---

## Acknowledgments

This work was performed using computational resources at the UNC Longleaf High Performance Computing Cluster. We thank the Phanstiel Lab and collaborators for feedback and support.

**Funding:** [Funding sources TBD]

---

## Contact

**JP Flores**  
PhD Candidate, Bioinformatics & Computational Biology  
University of North Carolina at Chapel Hill  
Email: jpflores@unc.edu

**Lab:** Phanstiel Lab (https://phanstiel-lab.med.unc.edu/)

For questions about:
- **Code/Analysis:** Open a GitHub issue or email JP Flores
- **Experimental methods:** Contact [Lab PI TBD]
- **Data access:** See GEO accession GSE_TBD or contact JP Flores

---

## Version History

See `git log` for detailed commit history.

**Major releases:**
- v1.0.0 (TBD): Initial manuscript submission
- v1.1.0 (TBD): Post-review revisions
- v2.0.0 (TBD): Final published version

---

**Last updated:** November 2024
