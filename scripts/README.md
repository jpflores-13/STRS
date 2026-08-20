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
| `rad21Peaks_lost_gained_homerInput.R` | RAD21 narrowPeak (control, sorbitol) | HOMER input files | Split RAD21 peaks into retained (overlap sorbitol) vs. lost sets for motif analysis |
| `top500peaks_cutntag_homerInput.R` | CUT&Tag peaks | HOMER input files | Prepare top 500 peaks for each protein for motif analysis |
| `atacPeaks_lost_gained_homerInput.R` | ATAC-seq peaks, differential loops | HOMER input files | Identify ATAC peaks overlapping gained/lost loop anchors, export BED for motif analysis |
| `extract_gained_anchor_promoters_homerInput.R` | Differential loops, TxDb hg38 | HOMER input files, gene list | Identify gene promoters (TSS ± 2000/500 bp) overlapping gained loop anchors |
| `gained_promoter_motif_homerInput_SP1.R` | Differential loops, RNA-seq LRT results | HOMER input files, motif enrichment plot | Test whether gained-anchor promoters are enriched for SP/KLF motifs vs. an expression-matched background (`nullranges::matchRanges()`) |

### EISA Processing (Exon-Intron Split Analysis)

| Script | Input | Output | Description |
|--------|-------|--------|-------------|
| `eisa_bam_processing.R` | HISAT2-aligned BAM files (RNA-seq timecourse) | `exon_counts.rds`, `intron_counts.rds` | EISA Part 1: count exonic/intronic reads via `Rsubread::featureCounts` (16 threads, ~30–40 min) |
| `eisa_auxinRAD21_bam_processing.R` | Per-replicate HISAT2-aligned BAMs (HCT116 mAID2-RAD21, n=6) | `EISA/per_replicate/exon_counts.rds`, `intron_counts.rds` | EISA Part 1 for the RAD21 auxin-degron dataset; per-replicate counts so eisaR/edgeR can estimate dispersion with 2 reps/condition |
| `eisa_auxinRAD21_bam_processing.sh` | — | SLURM job | Longleaf submission wrapper for `eisa_auxinRAD21_bam_processing.R` (`sbatch eisa_auxinRAD21_bam_processing.sh`) |

### Additional Processing

| Script | Input | Output | Description |
|--------|-------|--------|-------------|
| `hic_timecourse_loop_enrichment.R` | Timecourse Hi-C, loops | Enrichment data | Calculate loop enrichment across timecourse |

### External Data

| Script | Input | Output | Description |
|--------|-------|--------|-------------|
| `encode_chipseq_download.sh` | `data/external/encode_chipseq/metadata.tsv` (ENCODE batch download) | `metadata_filtered.tsv`, `beds/*.bed.gz` | SLURM job; filters ENCODE metadata to optimal IDR-thresholded (TF) and replicated (histone) GRCh38 peaks for HEK293/HEK293T, downloads passing BED files |

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
| `surveyPlot_gainedLoops_degron_10kb.R` | Gained loops in CTCF/RAD21 degron backgrounds |
| `surveyPlot_gainedLoops_RAD21_degron_10kb.R` | Gained loops with RAD21 degron Hi-C |
| `surveyPlot_gainedLoops_HCT116_WT_10kb.R` | Gained loops in HCT116 wild-type cells |
| `surveyPlot_gainedLoopGenes_rnaTimecourse_10kb.R` | Gained loop genes with RNA-seq timecourse overlay |
| `surveyPlots_24h_peaking.R` | Survey plots for loops peaking at 24 h |
| `surveyPlot_HiC_CutnTag.R` | Hi-C with CUT&Tag tracks |
| `surveyPlot_HiC_CutnTag_CTCF_RAD21.R` | Hi-C with CTCF and RAD21 tracks |
| `surveyPlot_compartments.R` | Compartment eigenvector tracks at example loci |
| `surveyPlot_compartments_YAPdTAD.R` | Compartment tracks at YAP1-dependent dTAD locus |

### ATAC-seq Analyses

| Script | Description |
|--------|-------------|
| `atac_clustering_apa.R` | APA analysis of Hi-C signal at ATAC peak clusters |
| `homer_atacAnchors.sh` | Run HOMER motif enrichment on ATAC peaks at loop anchors |
| `MotifEnrichment_CTCF_RAD21_ATAC.R` | Compare motif enrichment across CTCF, RAD21, and ATAC peak sets |

### EISA (Exon-Intron Split Analysis)

| Script | Description |
|--------|-------------|
| `eisa_analysis.R` | EISA of RNA-seq timecourse to separate transcriptional from post-transcriptional changes |
| `eisa_lineplots.R` | Line plots of EISA scores across timepoints |
| `eisa_scatter.R` | Scatter plots of exonic vs. intronic log2FC |
| `eisa_auxinRAD21_analysis.R` | EISA for RAD21 auxin degron RNA-seq data |
| `eisa_auxinRAD21_barplot.R` | Barplots of EISA results for RAD21 degron |
| `eisa_auxinRAD21_histogram.R` | Histograms of EISA scores for RAD21 degron |

### Domain and Compartment Analyses

| Script | Description |
|--------|-------------|
| `domain_analysis.sh` | Compute TAD insulation scores from Hi-C contact maps |
| `domain_analysis_visualizations.R` | Visualize insulation score and domain boundary changes |
| `domainAnalysis_YAPdTAD.R` | TAD analysis at YAP1-dependent domain loci |
| `domain_analysis_YAPdTAD.sh` | Compute insulation at YAP1-dependent dTAD regions |
| `compartmentStrength_boxplot.R` | Boxplots of compartment strength (saddle scores) across conditions |
| `compartment_strength_perRep.sh` | Per-replicate compartment strength calculation |

### Auxin Degron Analyses

| Script | Description |
|--------|-------------|
| `auxinDegronQuant.R` | Quantify depletion efficiency from CellProfiler output |
| `auxinRAD21_gainedLoops_barplot.R` | Barplots of gained loop counts in RAD21 degron |
| `auxinRAD21_gainedLoops_histogram.R` | Histograms of loop scores in RAD21 degron vs. control |
| `auxinRAD21_locusPlot.R` | Example locus plots for RAD21 degron Hi-C |
| `run_cellprofiler.sh` | Shell wrapper to run CellProfiler on auxin degron microscopy images |
| `CTCF_RAD21_auxin_degron.cppipe` | CellProfiler pipeline for CTCF/RAD21 degron image analysis |

### Additional Analyses

| Script | Description |
|--------|-------------|
| `anchor_expression_summary.R` | Summarize gene expression at loop anchor regions |
| `active_promoter_enrichment.R` | Enrichment of active promoter marks at loop anchors |
| `sorbitol_genes_at_loop_anchors.R` | Identify and characterize sorbitol-responsive genes at anchors |
| `temporal_trajectory_lineplot.R` | Line plots of expression trajectories for anchor genes |
| `temporal_trajectory_lineplot_allExpressedGenes.R` | Trajectory plots for all expressed genes (background comparison) |
| `sox8_counts.R` | Count-level analysis for SOX8 locus as an example |
| `encode_overlap_scatter.R` | Scatter of loop anchor overlap with ENCODE ChIP-seq datasets |
| `homerQQ.R` | QQ scatter plots of HOMER motif enrichment (-log10 p-value) |
| `homer_gained_promoters.sh` | HOMER motif analysis on gained-loop promoter anchors |
| `homer_gained_promoters_SP1.sh` | HOMER analysis focused on SP1/KLF motifs at gained promoters |
| `homer_rad21Anchors.sh` | HOMER motif analysis on RAD21 CUT&Tag peaks at loop anchors |
| `visualize_gained_promoter_motifs.R` | Visualize HOMER motif enrichment at gained-loop promoters |
| `DifferentialMatrix.R` | Compute differential contact matrices between conditions |
| `ctcf_retentionFigure.R` | Composite figure for CTCF retention statistics |
| `cre_cre_apa.R` | APA analysis between pairs of CRE peaks (CRE-CRE loops) |

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
| `FigureS2.R` | `figures/FigureS2.pdf` | Compartment and TAD domain organization under hyperosmotic stress |
| `FigureS4.R` | `figures/FigureS4.pdf` | Extended CUT&Tag and RNA-seq analyses |
| `FigureS5.R` | `figures/FigureS5.pdf` | Additional Hi-C validation and controls |
| `FigureS6.R` | `figures/FigureS6.pdf` | CUT&Tag binding scatter plots for CTCF, RAD21, and YAP1 (control vs. sorbitol) |
| `FigureS8.R` | `figures/FigureS8.pdf` | CTCF retention: log2FC vs. CPM bins, PWM score, RAD21 FC, and promoter density |
| `FigureS9.R` | `figures/FigureS9.pdf` | H3K27ac differential analysis: MA plot, anchor overlap bar plot, and anchor vs. between-anchor density |
| `FigureS10.R` | `figures/FigureS10.pdf` | CTCF and RAD21 auxin degron validation (microscopy images and GFP quantification) |
| `FigureS11.R` | `figures/FigureS11.pdf` | Enhancer-promoter loop enrichment and CRE-CRE APA for H3K27ac peaks |
| `FigureS12.R` | `figures/FigureS12.pdf` | HOMER motif QQ plots and ENCODE ChIP-seq overlap for CTCF, RAD21, and ATAC anchors |
| `FigureS13.R` | `figures/FigureS13.pdf` | DoG transcription at select loci (modest effect); RNA-seq signal tracks for SAMD4A and KDM6A |

> **Note:** `FigureS3.pdf` and `FigureS7.pdf` are collaborator-provided figures (nuclear volume quantification and CTCF/RAD21 fluorescence redistribution imaging, respectively); there is no corresponding R script for either.

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
| `calculate_apa_score.R` | Computes APA score (median foreground / median background pixel ratio) from a normalized APA matrix |
| `ggplot2_pgTheme.R` | Custom ggplot2 theme for consistent figure styling |
| `make_norm_matrix.R` | Matrix normalization functions for APA analysis |
| `rnaseqTimecourse_helper_functions.R` | Helper functions for RNA-seq timecourse GO enrichment (e.g., `compare_go_across_timepoints()`) |
| `validate_plotgardener_genes.R` | Gene annotation validation for plotgardener |
| `setup_cooltools_env.sh` | One-time interactive setup of a conda env (`cooltools_env`) with `hic2cool` and `cooltools`; not sourced by R scripts — run manually on a Longleaf login node |

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
