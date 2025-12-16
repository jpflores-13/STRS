#!/bin/bash
#SBATCH --job-name=homer_gained_promoters
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_80
#SBATCH --mail-user=jflores@unc.edu
#SBATCH --output=lslogs/homer_gained_promoters_%j.out
#SBATCH --error=logs/homer_gained_promoters_%j.err

# Load HOMER
ml homer/5.1

# Set working directory
cd /work/users/j/p/jpflores/projects/STRS

# Create output directories
mkdir -p data/processed/cutntag/homer_motifs/gained_anchor_promoters
mkdir -p data/processed/cutntag/homer_preparsed
mkdir -p logs

echo "========================================="
echo "HOMER Motif Analysis: Gained Anchor Promoters"
echo "========================================="
echo ""

# Print input file stats
echo "Input files:"
echo "  Foreground: $(wc -l < data/processed/cutntag/homer_input/gained_anchor_promoters.bed) promoters"
echo "  Background: $(wc -l < data/processed/cutntag/homer_input/all_promoters.bed) promoters"
echo ""

# Run HOMER with all promoters as background
echo "Analyzing promoters at gained loop anchors vs all promoters..."
findMotifsGenome.pl \
  data/processed/cutntag/homer_input/gained_anchor_promoters.bed \
  hg38 \
  data/processed/cutntag/homer_motifs/gained_anchor_promoters \
  -size given \
  -bg data/processed/cutntag/homer_input/all_promoters.bed \
  -mask \
  -p 8 \
  -preparsedDir data/processed/cutntag/homer_preparsed

echo ""
echo "========================================="
echo "Analysis complete!"
echo "========================================="
echo "Results saved to:"
echo "  data/processed/cutntag/homer_motifs/gained_anchor_promoters/homerResults.html"