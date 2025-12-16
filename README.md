# Hyperosmotic stress induces a large-scale rewiring of 3D chromatin interactions 

This repository contains code, analysis scripts, and figure generation for our study investigating the role of 3D chromatin structure in response to hyperosmotic stress.

> **Flores JP**, Perreault AA, Drum Z, Xu C, Cruz Alonso D, Petros G, Wu Y, Quiroga-Barber IY, Sahasrabudhe I, Demmerle J, Wang GG, Cai D, Phanstiel DH  
> *Hyperosmotic stress induces a large-scale rewiring of 3D chromatin interactions *  
> Preprint: **TBD DOI** · Journal: **TBD**

## Overview

We integrated Hi-C, CUT&Tag, and RNA-seq data to characterize how sorbitol treatment (hyperosmotic stress) affects chromatin loop formation and gene expression in HEK293 and HCT116 cells. Our analysis reveals large-scale rewiring of chromatin interactions with distinct patterns of loop gain and loss, driven by differential CTCF and cohesin (RAD21) retention at loop anchors.

**HIGHLIGHTS:**
- Acute hyperosmotic stress causes a global  loss of existing chromatin loops accompanied by the formation of hundreds of de novo loops.
- Gained loop anchors are enriched for retained CTCF and cohesin.
- CTCF retention is favored at high-occupancy, cohesin-stabilized, promoter-proximal sites with strong CTCF motifs.
- Loop remodeling is largely decoupled from genome-wide transcription but can associate with both activation and repression at specific loci.

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
- Peak calling: MACS2
- Differential analysis: DESeq2
- Proteins analyzed: CTCF, RAD21, YAP1, H3K27ac
- Motif analysis: HOMER

### RNA-seq Analysis
- Alignment: HiSat2
- Quantification: featureCounts
- Differential expression: DESeq2 with likelihood ratio test (LRT)
- Clustering: K-means on variance-stabilized counts

### Software Versions
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
- **Figure 1:** Hyperosmotic stress induces complete rewiring of chromatin interactions.
- **Figure 2:** Sorbitol-induced loops are more punctate, form weaker chromatin domains, and peak at 1 hour of treatment.
- **Figure 3:** CTCF and cohesin are retained at gained loop anchors following hyperosmotic stress.
- **Figure 4:** High-occupancy, strong-motif, cohesin-stabilized, promoter-proximal sites preferentially retain CTCF during hyperosmotic stress.
- **Figure 5:** Sorbitol-induced loops are not strongly associated with transcriptional changes.

### Supplementary Figures
- **Figure S1:** Hyperosmotic stress induces predominantly de novo looping that frequently reuses pre-existing anchor sites and favors long-range interactions.
- **Figure S2:** Covariate matching controls for loop size and interaction frequency in gained loop comparisons.
- **Figure S3:** CTCF and cohesin binding is selectively retained at anchors of sorbitol-induced chromatin loops.

---

## Citation

If you use code or data from this repository, please cite:

```
Flores JP, et al. (2024)
Hyperosmotic stress induces chromatin loop reorganization through CTCF retention
[Journal TBD] [DOI TBD]
```

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

## Acknowledgements

We thank Erika Deoudes for data visualization, illustration, proofreading, and typesetting. We thank Samantha Pattenden for use of the Covaris LE220 instrument, which was provided by the North Carolina Biotechnology Center Institutional Development Program grant 2017-IDG-1005. We also thank Brian Golitz and the UNC CRISPR Core for technical assistance.

---

---

## Funding

This work was supported in part by the Howard Hughes Medical Institute (Gilliam Fellows Program #GT16825 to J.P.F.), the National Institutes of Health (R35GM128645 to D.H.P.; R01CA271603 to D.H.P. and G.W.), and the Department of Defense Kidney Cancer Idea Development Award (W81XWH2210900 to D.C.). Z.A.D. was supported by the Seeding Postdoctoral Innovators in Research and Education (SPIRE) Postdoctoral Training Program. D.C.A. and G.P. were supported by the Postbaccalaureate Research Education Program (PREP). D.C. was supported by the Department of Defense Kidney Cancer Idea Development Award (W81XWH2210900, D.C.) and the National Institutes of Health (R35GM142837, D.C.). J.D. was supported by the National Cancer Institute (NCI) training grant T32CA009110. A.A.P. was supported by the Cancer Epigenetics Training Program (5T32-CA217824) and an Elon University Faculty Research & Development grant. I.Y.Q.-B. was supported by a BrightFocus Foundation Fellowship (Fellowship 911831). 

---

## Contact

**JP Flores**  
PhD Candidate, Bioinformatics & Computational Biology  
University of North Carolina at Chapel Hill  
Email: jflores@unc.edu

**Lab:** Phanstiel Lab (https://phanstiel-lab.med.unc.edu/)

For questions about:
- **Code/Analysis:** Open a GitHub issue or email JP
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

**Last updated:** November 11, 2025
