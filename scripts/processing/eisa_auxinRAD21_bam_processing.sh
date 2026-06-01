#!/bin/bash

# ##############################################################################
# filename:    eisa_auxinRAD21_bam_processing.sh
# author:      JP Flores
# date:        2026-05-29
# project:     STRS
# description: SLURM submission script for eisa_auxinRAD21_bam_processing.R.
#              Runs featureCounts on the 6 per-replicate AuxinRAD21 BAM files
#              to generate exonic and intronic count matrices for eisaR EISA.
#              Runtime ~30-40 minutes. Requests 16 CPUs and 64GB RAM to match
#              the featureCounts nthreads = 16 setting.
# usage:       sbatch eisa_auxinRAD21_bam_processing.sh
# ##############################################################################

#SBATCH --job-name=eisa_auxinRAD21_bam
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --output=/work/users/j/p/jpflores/projects/STRS/logs/eisa_auxinRAD21_bam_%j.out
#SBATCH --error=/work/users/j/p/jpflores/projects/STRS/logs/eisa_auxinRAD21_bam_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jflores@unc.edu

## ── Environment --------------------------------------------------------------

module purge
module load r/4.5.0
module load samtools

## ── Run ----------------------------------------------------------------------

cd /work/users/j/p/jpflores/projects/STRS

echo "Job started: $(date)"
echo "Running on node: $(hostname)"
echo "SLURM job ID: ${SLURM_JOB_ID}"

Rscript scripts/processing/eisa_auxinRAD21_bam_processing.R

echo "Job finished: $(date)"