#!/bin/bash

# encode_chipseq_download.sh -----------------------------------------------
# Author:      JP Flores
# Date:        2026-04-16
# Project:     STRS
# Description: Filters metadata.tsv from ENCODE batch download to retain
#              only optimal IDR thresholded peaks (TF ChIP-seq) and
#              replicated peaks (Histone ChIP-seq) in GRCh38 for HEK293
#              and HEK293T cell lines. Downloads all passing BED files
#              into a dedicated beds/ subdirectory.
# Input:       data/external/encode_chipseq/metadata.tsv
# Output:      data/external/encode_chipseq/metadata_filtered.tsv
#              data/external/encode_chipseq/beds/*.bed.gz
# --------------------------------------------------------------------------

#SBATCH --job-name=encode_download
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jflores@unc.edu
#SBATCH --output=logs/encode_chipseq_download_%j.out
#SBATCH --error=logs/encode_chipseq_download_%j.err

# --------------------------------------------------------------------------
# Parameters
# --------------------------------------------------------------------------

PROJECT_DIR="/work/users/j/p/jpflores/projects/STRS"
ENCODE_DIR="${PROJECT_DIR}/data/external/encode_chipseq"
METADATA="${ENCODE_DIR}/metadata.tsv"
METADATA_FILTERED="${ENCODE_DIR}/metadata_filtered.tsv"
BEDS_DIR="${ENCODE_DIR}/beds"

# --------------------------------------------------------------------------
# Setup
# --------------------------------------------------------------------------

mkdir -p "${BEDS_DIR}"
mkdir -p "${PROJECT_DIR}/logs"

echo "=========================================="
echo "ENCODE ChIP-seq Download"
echo "Started: $(date)"
echo "=========================================="

# --------------------------------------------------------------------------
# Filter metadata.tsv
#
# Keep rows where:
#   TF ChIP-seq  : col2 == "bed narrowPeak"
#                  col5 == "optimal IDR thresholded peaks"
#                  col6 == "GRCh38"
#                  col8 == "TF ChIP-seq"
#
#   Histone ChIP : col2 == "bed narrowPeak"
#                  col5 == "replicated peaks"
#                  col6 == "GRCh38"
#                  col8 == "Histone ChIP-seq"
# --------------------------------------------------------------------------

echo "Filtering metadata.tsv..."

## Write header first
head -1 "${METADATA}" > "${METADATA_FILTERED}"

## Append TF ChIP-seq rows (optimal IDR thresholded peaks)
awk -F'\t' '
  NR > 1 &&
  $2 == "bed narrowPeak" &&
  $5 == "optimal IDR thresholded peaks" &&
  $6 == "GRCh38" &&
  $8 == "TF ChIP-seq"
' "${METADATA}" >> "${METADATA_FILTERED}"

## Append Histone ChIP-seq rows (replicated peaks)
awk -F'\t' '
  NR > 1 &&
  $2 == "bed narrowPeak" &&
  $5 == "replicated peaks" &&
  $6 == "GRCh38" &&
  $8 == "Histone ChIP-seq"
' "${METADATA}" >> "${METADATA_FILTERED}"

## Report how many files passed the filter
N_FILTERED=$(awk 'NR > 1' "${METADATA_FILTERED}" | wc -l)
echo "Files passing filter: ${N_FILTERED}"
echo "Filtered metadata saved to: ${METADATA_FILTERED}"

# --------------------------------------------------------------------------
# Download BED files
#
# Extract column 48 (File download URL) from filtered metadata,
# skip the header, and download each file into beds/
# --------------------------------------------------------------------------

echo ""
echo "Downloading ${N_FILTERED} BED files to ${BEDS_DIR}..."

awk -F'\t' 'NR > 1 {print $48}' "${METADATA_FILTERED}" | \
  xargs -n 1 -P 4 curl -L --silent --output-dir "${BEDS_DIR}" -O

echo ""
echo "Download complete: $(date)"
echo "Files in beds/: $(ls "${BEDS_DIR}" | wc -l)"
echo "=========================================="