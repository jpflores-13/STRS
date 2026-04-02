#!/bin/bash
#SBATCH --job-name=homer_atac_gained_lost
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=homer_atac_gained_lost_%j.out
#SBATCH --error=homer_atac_gained_lost_%j.err

# Load HOMER
module load homer/5.1

# Set working directory
cd /work/users/j/p/jpflores/projects/STRS

# Create output directories
mkdir -p data/processed/atac/homer_motifs/atac_gained_anchors
mkdir -p data/processed/atac/homer_motifs/atac_lost_anchors
mkdir -p data/processed/atac/homer_preparsed

echo "========================================="
echo "HOMER Motif Analysis: ATAC Peaks at Gained vs Lost Loop Anchors"
echo "========================================="
echo ""

# Run HOMER on ATAC peaks at gained anchors
echo "Analyzing ATAC peaks at GAINED loop anchors..."
findMotifsGenome.pl \
  data/processed/atac/homer_input/gained_anchor_atac.bed \
  hg38 \
  data/processed/atac/homer_motifs/atac_gained_anchors \
  -size 200 \
  -mask \
  -p 8 \
  -preparsedDir data/processed/atac/homer_preparsed

echo "Gained anchor ATAC peaks analysis complete!"
echo ""

# Run HOMER on ATAC peaks at lost anchors
echo "Analyzing ATAC peaks at LOST loop anchors..."
findMotifsGenome.pl \
  data/processed/atac/homer_input/lost_anchor_atac.bed \
  hg38 \
  data/processed/atac/homer_motifs/atac_lost_anchors \
  -size 200 \
  -mask \
  -p 8 \
  -preparsedDir data/processed/atac/homer_preparsed

echo "Lost anchor ATAC peaks analysis complete!"
echo ""

# Summary
echo "========================================="
echo "Summary"
echo "========================================="
echo "Results saved to:"
echo "  Gained: data/processed/atac/homer_motifs/atac_gained_anchors/homerResults.html"
echo "  Lost: data/processed/atac/homer_motifs/atac_lost_anchors/homerResults.html"

echo ""
echo "Job complete!"