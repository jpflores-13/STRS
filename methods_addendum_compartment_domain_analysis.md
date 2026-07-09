# Methods addendum: compartment and chromatin domain analysis

This note documents terminology and methods for the A/B compartment and
chromatin domain (insulation) analyses added in `scripts/analysis/domain_analysis.sh`,
`domain_analysis_YAPdTAD.sh`, `domainAnalysis_YAPdTAD.R`,
`domain_analysis_visualizations.R`, `compartmentStrength_boxplot.R`, and
`compartment_strength_perRep.sh` (feeding `FigureS2.R` / `FigureS2.pdf`). This
analysis is not yet described in `v1_STRS_Paper.pdf` — the text below is drafted
in the same style as the existing Methods section and is intended to be pasted
into the manuscript's Methods (under a new "Compartment and chromatin domain
analysis" subsection, after "Aggregate Hi-C analysis") once the source
manuscript document is updated.

## Terminology (for consistency with the paper and the grant)

- **"TAD" is already used in this manuscript** for YAP1's intrinsically
  disordered transcriptional/transcription activation domain (the
  eGFP-YAP1ΔTAD dominant-negative construct, e.g. paper.txt lines 623-625,
  656, 1050). To avoid collision, the manuscript never uses "TAD" for the
  topologically-associating-domain-like structures called from Hi-C
  insulation scores. Those are instead called **"chromatin domains"**
  ("domain boundaries," "domain-centric" when referring to the block of
  contacts inside a loop) — see paper.txt lines 77, 497, 501, 525, 1311.
- **A/B compartments** are called genome compartments / "A compartment" /
  "B compartment," derived from the sign of the leading eigenvector (PC1) of
  the Hi-C contact matrix. This terminology is standard and does not collide
  with anything else in the paper.
- Several script/variable names still use "TAD" internally
  (`domain_analysis.sh`, `aggreTAD.R`, `aggregateTAD.R`) — these predate the
  naming decision above and are code-internal only; the manuscript-facing
  text should always use "chromatin domain," never "TAD," for these
  structures.

## Draft Methods text

> **Compartment and chromatin domain analysis**
>
> A/B compartments were called from ICE-balanced Hi-C contact matrices at
> 10-kb resolution using the leading eigenvector (PC1) of the intrachromosomal
> correlation matrix (`cooltools eigs-cis`). PC1 sign was oriented per
> chromosome using a GC-content phasing track such that positive values
> correspond to the A (active, GC-rich) compartment and negative values to
> the B compartment. Compartment strength was quantified from saddle plots
> (`cooltools saddle`, 50 bins) generated from the control-sample eigenvector
> and condition-matched expected contact frequency, and summarized per
> replicate as (AA + BB) / (AB + BA) using the corner bins of the saddle
> matrix.
>
> Chromatin domain boundaries were identified from ICE-balanced 10-kb
> matrices using a diamond insulation score with a 250-kb sliding window
> (`cooltools insulation`). Mean insulation score profiles were computed
> across called domain boundaries to compare boundary strength between
> control and sorbitol-treated conditions [and, separately, between
> eGFP-YAP1 and eGFP-YAP1ΔTAD genotypes].
>
> *(Tool citation needed: Abdennur & Mirny, cooltools, Bioinformatics 2020 /
> Open2C et al. 2024, in addition to the existing Rao et al. 2014 Hi-C
> citation already used for library generation.)*

## Answer to "how did we call compartments and domains" (for the grant)

- **Compartments:** "A/B compartments," from Hi-C eigenvector (PC1)
  decomposition; "compartment strength" for the saddle-plot summary.
- **Domains:** "chromatin domains" / "domain boundaries" — deliberately
  **not** "TADs," because "TAD" in this project already refers to YAP1's
  transcriptional activation domain (eGFP-YAP1ΔTAD). Use "chromatin domain"
  or "insulation-defined domain" in the grant to stay consistent with the
  paper.
