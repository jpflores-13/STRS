# Scripts Directory

This directory contains all analysis code for the STRS project, organized into four main categories: processing raw data, performing downstream analyses, generating manuscript figures, and utility functions.

---

##  Directory Structure

```
scripts/
├── processing/    # Raw data → processed outputs (loops, APAs, peaks)
├── analysis/      # Downstream statistical analyses and visualizations
├── figures/       # Final manuscript figure generation
└── utils/         # Helper functions and themes
```

---

## Processing Scripts (`processing/`)

These scripts transform raw data into processed outputs for downstream analysis. They are **computationally intensive** and should be run on HPC compute nodes.

### Core Hi-C Processing

| Script | Input | Output | Description |
|--------|-------|--------|-------------|
| `extractLoops.R` | Loop calls (BEDPE), Hi-C maps | `sipLoops_*.rds` | Merge SIP loop calls and extract pixel counts from .hic files |
| `callDiffLoops.R` | Loop counts | `diffLoops_*.rds` | DESeq2 differential loop analysis (control vs sorbitol) |

### APA Matrix Calculations

| Script | Input | Output | Description |
|--------|-------|--------|-------------|
| `calcAPA_gained.R` | Differential loops, Hi-C maps | 8 APA matrices | Aggregate peak analysis for gained loops across 4 genotypes × 2 conditions |
| `calcAPA_lost.R` | Differential loops, Hi-C maps | 8 APA matrices | Aggregate peak analysis for lost loops across 4 genotypes × 2 conditions |
| `calcAPA_timecourse_gained.R` | Differential loops, timecourse Hi-C | 8 APA matrices | Timecourse APA for gained loops (0h-24h) |
| `calcAPA_timecourse_lost.R` | Differential loops, timecourse Hi-C | 8 APA matrices | Timecourse APA for lost loops (0h-24h) |

### HOMER Preparation

| Script | Input | Output | Description |
|--------|-------|--------|-------------|
| `ctcfPeaks_lost_gained_homerInput.R` | CTCF peaks, differential loops | HOMER input files | Prepare peak lists for motif analysis at gained/lost loop anchors |
| `top500peaks_cutntag_homerInput.R` | CUT&Tag peaks | HOMER input files | Prepare top 500 peaks for each protein for motif analysis |

### Additional Processing

| Script | Input | Output | Description |
|--------|-------|--------|-------------|
| `hic_timecourse_loop_enrichment.R` | Timecourse Hi-C, loops | Enrichment data | Calculate loop enrichment across timecourse |

---

##  Analysis Scripts (`analysis/`)

Downstream statistical analyses and exploratory visualizations. These produce data for figures and generate plots in `plots/` directory.

### Hi-C Loop Analysis

| Script | Description |
|--------|-------------|
| `diffLoops_MAplot.R` | MA plot for differential loop analysis |
| `loop_size_comparison.R` | Compare sizes of gained vs lost loops |
| `loop-type-enrichment.R` | Enrichment analysis for different loop categories |
| `strengthened_weakened_loop_analysis.R` | Analysis of pre-existing loops that strengthen or weaken |
| `exampleRegion_gained.R` | Generate example genomic regions with gained loops |
| `exampleRegion_lost.R` | Generate example genomic regions with lost loops |
| `exampleRegion_gainedLost.R` | Generate regions showing both gained and lost loops |
| `normalizedAPAs_genotype.R` | Compare normalized APAs across genotypes / cell types |
| `aggreTAD.R` | Aggregate TAD analysis |

### Hi-C Timecourse

| Script | Description |
|--------|-------------|
| `hic_timecourse.R` | General timecourse Hi-C analysis |
| `hic_timecourse_apa_enrichment.R` | APA enrichment over timecourse |
| `hic_timecourse_loop_enrichment_line_plot.R` | Line plots of APA score for loops |

### CUT&Tag Analyses

| Script | Description |
|--------|-------------|
| `cutntag_differentialAnalysis.R` | Differential binding analysis for CUT&Tag data |
| `cutntag_MAplots.R` | MA plots for CUT&Tag differential analysis |
| `cutntag_HiCanchors_proportions.R` | Proportion of loop anchors with protein binding |
| `cutntag_HiCanchors_boxplots_signif.R` | Statistical testing of protein enrichment at anchors |
| `cutntag_HiCanchors_indiv_proteins.R` | Individual protein analysis at loop anchors |
| `cutntag_anchor_vs_between_boxplots.R` | Compare signal at anchors vs between loops |
| `cutntag_all_loops_anchor_vs_between_density.R` | Density plots for all loops |
| `cutntag_gained_loops_anchor_vs_between_density.R` | Density plots for gained loops specifically |
| `log2FC_vs_CPMbins_boxplots.R` | CPM binning strategy to avoid power bias for CTCF |

### CTCF-Specific Analyses

| Script | Description |
|--------|-------------|
| `ctcf_anchors_motifs.R` | CTCF motif analysis at loop anchors |
| `ctcf_binding_strength_control.R` | Control analysis for CTCF binding strength |
| `ctcf_retention_by_pwm_score.R` | CTCF retention as function of motif PWM score |
| `ctcf_retention_by_rad21_fc.R` | CTCF retention as function of RAD21 change |
| `ctcf_log2fc_comparison.R` | Compare CTCF log2FC across conditions |
| `ctcf_cooccupancy_pairedBarplots.R` | CTCF co-occupancy with other proteins |
| `CTCF_RAD21_grouped_barplot.R` | Grouped comparison of CTCF and RAD21 |
| `promoter_enrichment_ctcf_anchors.R` | Promoter enrichment at CTCF-bound anchors |
| `promoterEnrichment_ctcf_retained_lost.R` | Promoter enrichment by CTCF retention status |

### RNA-seq Timecourse Analyses

| Script | Description |
|--------|-------------|
| `rnaseqTimecourse_LRT.R` | Likelihood ratio test for timecourse RNA-seq |
| `rnaseqTimecourse_clusters_heatmap.R` | K-means clustering and heatmap of temporal patterns |
| `rnaseqTimecourse_loopAnchors_boxplots.R` | Gene expression at loop anchors over time |
| `GOanalysis_timecourse.R` | Gene ontology enrichment for temporal clusters |

### Enrichment Analyses

| Script | Description |
|--------|-------------|
| `promoter-enrichment.R` | General promoter enrichment at loop anchors |
| `ep-loop-enrichment.R` | Enhancer-promoter loop enrichment analysis |

### Motif Analysis (HOMER)

| Script | Description |
|--------|-------------|
| `homer_ctcfAnchors.sh` | Run HOMER on CTCF peaks at loop anchors |
| `homer_top500.sh` | Run HOMER on top 500 peaks for each protein |
| `protein_comparison_motifs.R` | Compare motif enrichment across proteins |

### Survey Plots (Multi-Track Visualizations)

| Script | Description |
|--------|-------------|
| `surveyPlot_gainedLoops_10kb.R` | Multi-track view of gained loops (Hi-C + tracks) |
| `surveyPlot_lostLoops_10kb.R` | Multi-track view of lost loops |
| `surveyPlot_gainedLoops_rnaTimecourse_10kb.R` | Gained loops with RNA-seq timecourse |
| `surveyPlot_lostLoops_rnaTimecourse_10kb.R` | Lost loops with RNA-seq timecourse |
| `surveyPlot_HiC_CutnTag.R` | Hi-C with CUT&Tag tracks |
| `surveyPlot_HiC_CutnTag_CTCF_RAD21.R` | Hi-C with CTCF and RAD21 tracks |

---

## Figure Generation Scripts (`figures/`)

Generate final manuscript figures as publication-quality PDFs. These scripts read processed data and produce outputs in `figures/`.

| Script | Output | Description |
|--------|--------|-------------|
| `Figure1.R` | `figures/Figure1.pdf` | Differential loop analysis and APA validation across cell lines |
| `Figure2.R` | `figures/Figure2.pdf` | Hi-C timecourse analysis of loop dynamics |
| `Figure3.R` | `figures/Figure3.pdf` | CUT&Tag enrichment at loop anchors (CTCF, RAD21, YAP1, H3K27ac) |
| `Figure4.R` | `figures/Figure4.pdf` | CTCF motif analysis and retention patterns |
| `Figure5.R` | `figures/Figure5.pdf` | RNA-seq timecourse and gene expression at loop anchors |
| `FigureS1.R` | `figures/FigureS1.pdf` | Loop overlap analysis and size distributions |
| `FigureS2.R` | `figures/FigureS2.pdf` | Additional Hi-C validation and controls |
| `FigureS3.R` | `figures/FigureS3.pdf` | Extended CUT&Tag and RNA-seq analyses |

**Usage:**
```bash
# Generate individual figures
Rscript scripts/figures/Figure1.R

# Or use Makefile
make fig1
make all  # Generate all figures
```

**Dependencies:** All figure scripts require processed data outputs from `processing/` and/or `analysis/` scripts.

---

## 🛠️Utility Functions (`utils/`)

Reusable helper functions sourced by other scripts. These should not be run directly.

| Script | Description |
|--------|-------------|
| `aggregateLoops.R` | Helper functions for aggregating loop data |
| `aggregateTAD.R` | Helper functions for TAD aggregation |
| `plotAggTAD.R` | Plotting functions for aggregate TAD analysis |
| `ggplot2_pgTheme.R` | Custom ggplot2 theme for consistent figure styling |
| `make_norm_matrix.R` | Matrix normalization functions for APA analysis |
| `validate_plotgardener_genes.R` | Gene annotation validation for plotgardener |

**Usage in scripts:**
```r
source("scripts/utils/ggplot2_pgTheme.R")
source("scripts/utils/aggregateLoops.R")
```

---

##  Typical Workflow

### 1. Processing
```bash
# Extract and call differential loops
Rscript scripts/processing/extractLoops.R
Rscript scripts/processing/callDiffLoops.R

# Calculate APA matrices
Rscript scripts/processing/calcAPA_gained.R
Rscript scripts/processing/calcAPA_lost.R
Rscript scripts/processing/calcAPA_timecourse_gained.R
Rscript scripts/processing/calcAPA_timecourse_lost.R
```

### 2. Analysis (Exploratory)
```bash
# Run specific analyses as needed
Rscript scripts/analysis/cutntag_HiCanchors_proportions.R
Rscript scripts/analysis/rnaseqTimecourse_clusters_heatmap.R

# These create plots in plots/ directory
```

### 3. Figure Generation (Iterative)
```bash
# Generate manuscript figures
make fig1
make fig2
# ... etc

# Or all at once
make
```
---

## Dependencies

### R Packages (managed via renv)
- **Bioconductor:** DESeq2, InteractionSet, mariner, plotgardener, GenomicRanges
- **Tidyverse:** dplyr, tidyr, ggplot2, stringr, purrr
- **Visualization:** RColorBrewer, ComplexHeatmap, cowplot
- **Statistics:** rstatix, clusterProfiler

See `renv.lock` in project root for exact versions.

### External Tools
- **HOMER:** Motif analysis (used by `homer_*.sh` scripts)
- **R ≥ 4.3.0**

---

## Notes

### Script Naming Conventions
- `*_analysis.R` → Statistical analyses
- `*_plot.R` or `*_boxplots.R` → Visualization scripts
- `calc*.R` → Computational/calculation scripts
- `surveyPlot_*.R` → Multi-track genomic visualizations
- `Figure*.R` → Final manuscript figures

### Output Locations
- **Processing scripts** → `data/processed/`
- **Analysis scripts** → `plots/` (exploratory)
- **Figure scripts** → `figures/` (final PDFs)
