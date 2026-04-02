#!/bin/bash
#SBATCH --job-name=homer_gained_promoters_SP1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=homer_gained_promoters_%j.out
#SBATCH --error=homer_gained_promoters_%j.err

# ── HOMER Motif Analysis: Gained Loop Anchor Promoters ──────────────────────
#
# Compares:
#   Focus      -> promoters (TSS ± 2kb/500bp) overlapping gained loop anchors
#   Background -> expression-matched promoters NOT at gained anchors
#                 (selected via matchRanges in gained_promoter_motif_homerInput.R)
#
# The explicit -bg flag tells HOMER to use our matched background instead of
# generating its own random genomic background — critical for this test.

module load homer/5.1

# ── Paths ────────────────────────────────────────────────────────────────────

PROJ_DIR="/work/users/j/p/jpflores/projects/STRS"
cd "$PROJ_DIR"

FOCUS_BED="data/processed/cutntag/homer_input/gained_anchor_promoters_focus_SP1.bed"
BG_BED="data/processed/cutntag/homer_input/gained_anchor_promoters_background_SP1.bed"
OUT_DIR="data/processed/cutntag/homer_motifs/gained_anchor_promoters_SP1"
PREPARSED_DIR="data/processed/cutntag/homer_preparsed"

mkdir -p "$OUT_DIR"
mkdir -p "$PREPARSED_DIR"

# ── Sanity checks ────────────────────────────────────────────────────────────

echo "========================================="
echo "HOMER Motif Analysis: Gained Loop Promoters"
echo "========================================="
echo ""

if [ ! -f "$FOCUS_BED" ]; then
  echo "ERROR: Focus BED not found: $FOCUS_BED"
  echo "Run gained_promoter_motif_homerInput_SP1.R first."
  exit 1
fi

if [ ! -f "$BG_BED" ]; then
  echo "ERROR: Background BED not found: $BG_BED"
  echo "Run gained_promoter_motif_homerInput_SP1.R first."
  exit 1
fi

FOCUS_N=$(wc -l < "$FOCUS_BED")
BG_N=$(wc -l < "$BG_BED")

echo "Focus regions:     $FOCUS_N"
echo "Background regions: $BG_N"
echo ""

# ── Run HOMER ────────────────────────────────────────────────────────────────
#
# Key flags:
#   -bg          : explicit background BED (our expression-matched set)
#   -size given  : use the actual coordinates from the BED file as-is,
#                  rather than re-centering/resizing around a summit.
#                  Promoter windows are already the right size from R.
#   -mask        : mask repeat regions (same as other analyses)
#   -p           : threads
#   -preparsedDir: reuse genome parsing cache across jobs

echo "Running findMotifsGenome.pl..."
echo ""

findMotifsGenome.pl \
  "$FOCUS_BED" \
  hg38 \
  "$OUT_DIR" \
  -bg "$BG_BED" \
  -size given \
  -mask \
  -p 8 \
  -preparsedDir "$PREPARSED_DIR"

echo ""
echo "========================================="
echo "Done!"
echo "Results in: $OUT_DIR"
echo ""
echo "Key files to check:"
echo "  $OUT_DIR/knownResults.txt   <- enrichment vs JASPAR/HOMER database"
echo "  $OUT_DIR/homerResults.txt   <- de novo motifs"
echo "  $OUT_DIR/homerResults.html  <- visual summary"
echo ""
echo "SP/KLF motifs to look for in knownResults.txt:"
echo "  SP1, SP2, SP3, SP4 (Specificity Protein family)"
echo "  KLF1-17            (Krüppel-Like Factor family)"
echo "  Both are GC-box binding factors with similar GGGCGG/CCGCCC consensus"
echo "========================================="