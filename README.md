# Hyperosmotic stress induces rewiring of 3D chromatin interactions

[![DOI](https://zenodo.org/badge/1094695813.svg)](https://doi.org/10.5281/zenodo.17989905)

This repository contains code, analysis scripts, and figure generation for our study investigating the role of 3D chromatin structure in response to hyperosmotic stress.

> **Flores JP**, Perreault AA, Drum Z, Xu C, Cruz Alonso D, Rashid T, Burgess JD, Fox GC, Petros G, Wu Y, Quiroga-Barber IY, Kim H, Sahasrabudhe I, Lee J, Black E, Li Y, Demmerle J, Strahl BD, Dowen JM, Wang GG, Cai D, Phanstiel DH  
> *Hyperosmotic stress induces rewiring of 3D chromatin interactions*  
> Preprint: https://www.biorxiv.org/content/10.64898/2025.12.19.695003v2 · Status: **Resubmitted; under review at *Cell Reports***

## Abstract

Cells continually encounter environmental stressors that challenge homeostasis. How three-dimensional (3D) chromatin structure contributes to these stress responses, particularly under hyperosmotic conditions, remains poorly understood. Here, using time-resolved Hi-C, CUT&Tag, auxin-inducible depletion, and RNA-seq, we map 3D chromatin structure, its molecular drivers, and transcriptional outcomes during the hyperosmotic stress response. Within 1 hour of sorbitol treatment, pre-existing loops and domains undergo genome-wide collapse, accompanied by the emergence of several hundred de novo, sorbitol-induced loops that are more punctate, longer-range, and transient. These newly formed loops weaken over time and largely dissipate by 24 hours, coincident with recovery of pre-existing chromatin structure. Loop reorganization is consistent across human cell types and hyperosmotic stimuli. CUT&Tag and degron experiments reveal that sorbitol-induced loops require cohesin but not CTCF. Newly formed loop anchors are enriched at active promoters containing SP and KLF family motifs. Genes located at these anchors show little immediate transcriptional change but are activated several hours after loop formation, consistent with loops functioning upstream of gene activation. Together, our findings show that hyperosmotic stress triggers a rapid, reversible, and CTCF-independent reorganization of 3D chromatin interactions that helps coordinate transcriptional adaptation.

**HIGHLIGHTS:**
- Hyperosmotic stress causes widespread loss and formation of chromatin loops.

- Stress-induced loop formation requires cohesin but not CTCF.

- Stress-induced loops are enriched at promoters with SP/KLF motifs.

- Genes at stress-induced loop anchors exhibit delayed expression

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

> **Note:** Large data files are not versioned in this repository.
> All raw and processed sequencing data generated in this study have been submitted to the NCBI Gene Expression Omnibus (GEO; https://www.ncbi.nlm.nih.gov/geo/). The Hi-C data are available under accession number GSE310051. The Hi-C data for HCT116-RAD21-mAID2 and HCT116-CTCF-mAID2 cells are available under accession number GSE312288. The RNA-seq (timecourse) data are available under accession number GSE310049. The CUT&Tag data are available under accession number GSE310047. The ATAC-seq data are available under accession number GSE329313. The RAD21-degron RNA-seq data are available under accession number GSE335237. See `data/*/README.md` for details.

---

## Quick Start

### Prerequisites

- R ≥ 4.3.0
- Required R packages (managed via `renv`)
- GNU Make (optional, for automated figure generation)

### Setup

```bash
# Clone the repository
git clone https://github.com/jpflores-13/STRS.git
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
make supplementary # Supplementary figures only (S1, S2, S4-S6, S8-S13; S3 and S7 have no scripts)
make fig1         # Generate Figure 1 only
make clean        # Remove generated figures
make help         # See all available commands
```

---

## Data Availability

### Raw Data
- **Hi-C:** FASTQ files in GEO: **GSE310051**
- **Degron Hi-C (HCT116-RAD21-mAID2 and HCT116-CTCF-mAID2):** FASTQ files in GEO: **GSE312288**
- **RNA-seq (timecourse):** FASTQ files in GEO: **GSE310049**
- **CUT&Tag:** FASTQ files in GEO: **GSE310047**
- **ATAC-seq:** FASTQ files in GEO: **GSE329313**
- **RAD21-degron RNA-seq:** FASTQ files in GEO: **GSE335237**

### Processed Data
All processed data are available in:
- Hi-C contact maps (`.hic` format)
- Loop calls (BEDPE format)
- Differential loop analysis results (`.rds`)
- Normalized APA matrices (`.rds`)
- CUT&Tag peak calls (narrowPeak format)
- CUT&Tag signal tracks (bigWig format)
- ATAC-seq peak calls (narrowPeak format)
- ATAC-seq signal tracks (bigWig format)
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

### ATAC-seq Analysis
- Alignment: BWA-MEM
- Peak calling: MACS2
- Differential analysis: DESeq2
- Used for chromatin accessibility context at loop anchors and motif enrichment

### RNA-seq Analysis
- Quantification: Salmon (v1.10.0) → tximeta/tximport → gene-level counts
- Differential expression: DESeq2, likelihood ratio test (LRT) for time-dependent genes;
  Wald tests + apeglm shrinkage for timepoint-specific contrasts
- Clustering: K-means on variance-stabilized, Z-scored counts

### Exon-Intron Split Analysis (EISA)
- Alignment: HISAT2 (v2.2.1)
- Quantification: featureCounts (Rsubread v2.24.0), exonic and intronic reads counted separately
- Purpose: distinguishes transcriptional (nascent, Δintron) from post-transcriptional
  (Δexon − Δintron) contributions to expression changes

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
- **Figure 1:** Hyperosmotic stress induces large-scale rewiring of chromatin interactions
- **Figure 2:** Sorbitol-induced loops are more punctate, form weaker chromatin domains, and peak at 1 hour of treatment
- **Figure 3:** CTCF, cohesin, and YAP1 are retained at gained loop anchors following hyperosmotic stress
- **Figure 4:** Sorbitol-induced chromatin loops require RAD21 and are enriched at promoters
- **Figure 5:** Sorbitol-induced loops are associated with transcriptional changes

### Supplementary Figures
- **Figure S1:** Hyperosmotic stress induces predominantly de novo looping that frequently reuses pre-existing anchor sites and favors long-range interactions
- **Figure S2:** 3D chromatin domain organization under hyperosmotic stress: compartments, saddle plots, and TAD insulation
- **Figure S3:** Hyperosmotic stress reduces nuclear volume (collaborator-provided figure; no R script in this repo)
- **Figure S4:** CTCF and cohesin binding is selectively retained at anchors of sorbitol-induced chromatin loops
- **Figure S5:** Covariate matching controls for loop size and interaction frequency in gained loop comparisons
- **Figure S6:** CUT&Tag binding scatter plots for CTCF, RAD21, and YAP1 in control vs. sorbitol conditions
- **Figure S7:** CTCF and RAD21 redistribute following hyperosmotic stress (collaborator-provided figure; no R script in this repo)
- **Figure S8:** CTCF retention analysis: log2FC vs. CPM bins, PWM score, RAD21 fold-change, and promoter density plots
- **Figure S9:** H3K27ac differential analysis at loop anchors: MA plot, anchor overlap bar plot, and anchor vs. between-anchor density
- **Figure S10:** CTCF and RAD21 auxin-inducible degron validation: microscopy images and nuclear GFP intensity quantification
- **Figure S11:** Enhancer-promoter loop enrichment and CRE-CRE APA analysis of H3K27ac peaks at gained loop anchors
- **Figure S12:** HOMER motif QQ scatter plots and ENCODE ChIP-seq overlap for CTCF, RAD21, and ATAC loop anchors
- **Figure S13:** Hyperosmotic stress induces downstream-of-gene (DoG) transcription at select loci

A graphical abstract submitted with the manuscript is available at `figures/GraphicalAbstract.pdf`.

---

## Citation

If you use code or data from this repository, please cite:

```
Flores JP, et al. (2026)
Hyperosmotic stress induces rewiring of 3D chromatin interactions
Resubmitted; under review at Cell Reports [DOI TBD]
```

For a specific version of this repository, cite the Zenodo DOI: doi.org/10.5281/zenodo.17989905

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

We thank Erika Deoudes for data visualization, illustration, proofreading, and typesetting. We thank Samantha Pattenden for use of the Covaris LE220 instrument, which was provided by the North Carolina Biotechnology Center Institutional Development Program grant 2017-IDG-1005. We also thank Brian Golitz and the UNC CRISPR Core for technical assistance. We thank Wendy Salmon and the UNC Hooker Imaging Core for assistance with fluorescence microscopy, including use of the GE IN Cell Analyzer 2200 and the Leica Stellaris 8 FALCON STED confocal microscope (NIH grant 1S10OD030300).

This work was supported in part by the Howard Hughes Medical Institute (Gilliam Fellows Program #GT16825 to J.P.F.), the National Institutes of Health (R35GM128645 to D.H.P.; R01CA271603 to D.H.P. and G.W.; R35GM142837 and R01CA303867 to D.C.; T32CA009110 to E.B.), and the Department of Defense Kidney Cancer Idea Development Award (W81XWH2210900 to D.C.). Confocal imaging at the UNC Hooker Imaging Core was supported by NIH grant 1S10OD030300. Z.A.D. was supported by the Seeding Postdoctoral Innovators in Research and Education (SPIRE) Postdoctoral Training Program at UNC Chapel Hill (5K12GM000678). J.D.B. was supported by a National Institute on Aging (NIA) K00 (K00 AG068509). T.R. was supported in part by a grant from the National Institute of General Medical Sciences (5T32 GM067553). G.C.F. was supported by a predoctoral fellowship from the National Human Genome Research Institute (F31HG14124-01). J.M.D. was supported by a grant from the National Institute of General Medical Sciences (R35GM152103). B.D.S. was supported by a grant from the National Institute of General Medical Sciences (R35GM126900). D.C.A. and G.P. were supported by the Postbaccalaureate Research Education Program (PREP) at UNC Chapel Hill (5R25GM089569). J.D. was supported by the National Cancer Institute (NCI) training grant T32CA009110. A.A.P. was supported by the Cancer Epigenetics Training Program (5T32-CA217824) and Elon University Faculty Research & Development grants. I.Y.Q.-B. was supported by a BrightFocus Foundation Fellowship (Fellowship 911831).

---

## Contact

**JP Flores, PhD**  
University of North Carolina at Chapel Hill  
Email: jflores@unc.edu (Alt: jpflores013@gmail.com)

**Lab:** Phanstiel Lab (https://phanstiel-lab.med.unc.edu/)

For questions about:
- **Code/Analysis:** Open a GitHub issue or email JP
- **Experimental methods:** Email JP
- **Data access:** The Hi-C data are available under accession number GSE310051. The Hi-C data for HCT116-RAD21-mAID2 and HCT116-CTCF-mAID2 cells are available under accession number GSE312288. The RNA-seq data are available under accession number GSE310049. The CUT&Tag data are available under accession number GSE310047. The ATAC-seq data are available under accession number GSE329313. The RAD21-degron RNA-seq data are available under accession number GSE335237. Contact JP Flores for other inquiries.

---

## Version History

See `git log` for detailed commit history.

**Major releases:**
- v1.0.0 (12/19/2025): Initial manuscript submission
- v2.0.0 (06/11/2026): Revised manuscript — updated title and author list, added ATAC-seq data (GSE329313), expanded supplementary figures (S6-S10, S12), added domain/compartment and auxin degron analyses
- v3.0.0 (07/30/2026): Resubmission to Cell Reports — added full manuscript PDF and graphical abstract, added Figure S3 (nuclear volume), fixed inverted Figure S12 description, merged Funding into Acknowledgements
- v3.1.0 (08/20/2026): Resubmission #3 — added Figure S7 (CTCF/RAD21 fluorescence redistribution imaging, collaborator-provided) per reviewer comment; renumbered Figures S7-S12 to S8-S13 accordingly (scripts, figure PDFs, Makefile, READMEs); updated Highlights and added eTOC blurb for Cell Reports submission

---

**Last updated:** August 20, 2026
