#!/bin/bash
#SBATCH --job-name=homer_ctcf_retained_lost
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=homer_ctcf_retained_lost_%j.out
#SBATCH --error=homer_ctcf_retained_lost_%j.err

# Load HOMER
module load homer/5.1

# Set working directory
cd /work/users/j/p/jpflores/projects/STRS

# Create output directories
mkdir -p data/processed/cutntag/homer_motifs/ctcf_retained_peaks
mkdir -p data/processed/cutntag/homer_motifs/ctcf_lost_peaks
mkdir -p data/processed/cutntag/homer_preparsed

echo "========================================="
echo "HOMER Motif Analysis: Retained vs Lost CTCF Peaks (Sorbitol Treatment)"
echo "========================================="
echo ""

# Run HOMER on retained CTCF peaks
echo "Analyzing RETAINED CTCF peaks (present in both control and sorbitol)..."
findMotifsGenome.pl \
  data/processed/cutntag/homer_input/ctcf_retained_peaks.bed \
  hg38 \
  data/processed/cutntag/homer_motifs/ctcf_retained_peaks \
  -size 200 \
  -mask \
  -p 8 \
  -preparsedDir data/processed/cutntag/homer_preparsed

echo "Retained peaks analysis complete!"
echo ""

# Run HOMER on lost CTCF peaks
echo "Analyzing LOST CTCF peaks (only in control, absent in sorbitol)..."
findMotifsGenome.pl \
  data/processed/cutntag/homer_input/ctcf_lost_peaks.bed \
  hg38 \
  data/processed/cutntag/homer_motifs/ctcf_lost_peaks \
  -size 200 \
  -mask \
  -p 8 \
  -preparsedDir data/processed/cutntag/homer_preparsed

echo "Lost peaks analysis complete!"
echo ""

# Summary
echo "========================================="
echo "Summary"
echo "========================================="
echo "Results saved to:"
echo "  Retained: data/processed/cutntag/homer_motifs/ctcf_retained_peaks/homerResults.html"
echo "  Lost: data/processed/cutntag/homer_motifs/ctcf_lost_peaks/homerResults.html"

echo ""
echo "Job complete!"
