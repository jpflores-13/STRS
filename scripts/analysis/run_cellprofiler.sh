#!/bin/bash
# =============================================================================
# run_cellprofiler.sh
# Author:      JP Flores
# Date:        2026-04-19
# Project:     STRS
# Description: Runs CellProfiler 4.2.8 headlessly on Longleaf using an
#              Apptainer container. This avoids compilation issues with
#              python-javabridge on RHEL9 + modern NumPy that prevent a
#              clean pip/conda install.
# Usage:       sbatch run_cellprofiler.sh
# =============================================================================
#
# ---- FIRST-TIME SETUP (run once interactively, not as a job) ---------------
#
#   module load apptainer
#   apptainer pull \
#       /work/users/j/p/jpflores/projects/STRS/containers/cellprofiler-4.2.8.sif \
#       docker://cellprofiler/cellprofiler:4.2.8
#
#   # Verify it works:
#   apptainer exec \
#       /work/users/j/p/jpflores/projects/STRS/containers/cellprofiler-4.2.8.sif \
#       python -m cellprofiler --version
#
# ----------------------------------------------------------------------------

#SBATCH --job-name=cellprofiler_auxin
#SBATCH --output=/work/users/j/p/jpflores/projects/STRS/scripts/analysis/logs/cellprofiler_%j.out
#SBATCH --error=/work/users/j/p/jpflores/projects/STRS/scripts/analysis/logs/cellprofiler_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16g
#SBATCH --time=2:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jflores@unc.edu

# --- Paths -------------------------------------------------------------------

PROJECT_DIR="/work/users/j/p/jpflores/projects/STRS"

## NOTE: the image directory name contains spaces (from IN Cell Analyzer output).
## A file list CSV is used instead of -i to ensure all 180 images are loaded.
## Generate with: find "IMAGE_DIR" -maxdepth 1 -name "*.tif" | sort | awk 'BEGIN{print "URI"} {print "file:" $0}' > image_file_list.csv
IMAGE_DIR="${PROJECT_DIR}/data/raw/imaging/auxin_images_flat"

PIPELINE="${PROJECT_DIR}/scripts/analysis/CTCF_RAD21_auxin_degron.cppipe"
OUTPUT_DIR="${PROJECT_DIR}/data/processed/imaging/cellprofiler_output"
LOG_DIR="${PROJECT_DIR}/scripts/analysis/logs"
SIF="${PROJECT_DIR}/containers/cellprofiler-4.2.8.sif"

# --- Environment -------------------------------------------------------------

module purge
module load apptainer

# --- Create output/log directories if they don't exist ----------------------

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"
mkdir -p "${PROJECT_DIR}/containers"

# --- Run CellProfiler headlessly via Apptainer ------------------------------
#
# apptainer exec: runs a command inside the container
# --bind: mounts your Longleaf paths into the container so it can see your
#         files. The container is isolated by default — without --bind, it
#         cannot access /work at all.
#
# CellProfiler flags:
#   -c             : run headless (no GUI)
#   -r             : run the pipeline
#   -p             : path to .cppipe file
#   --file-list    : explicit list of all image URIs (replaces -i, which only loads 1 image set headlessly)
#   -o             : output directory for CSVs
#
# CellProfiler is single-threaded by design — 1 CPU is correct.
# 16GB RAM is conservative but safe for 2048x2048 TIFFs.

apptainer exec \
    --bind "${PROJECT_DIR}:${PROJECT_DIR}" \
    "${SIF}" \
    python -m cellprofiler \
        -c \
        -r \
        -p "${PIPELINE}" \
        -o "${OUTPUT_DIR}"

echo "CellProfiler finished. Output written to: ${OUTPUT_DIR}"