#!/bin/bash

#SBATCH --job-name=domain_analysis
#SBATCH --output=/work/users/j/p/jpflores/projects/STRS/scripts/analysis/logs/domain_analysis_%j.out
#SBATCH --error=/work/users/j/p/jpflores/projects/STRS/scripts/analysis/logs/domain_analysis_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL,INCOMPLETE
#SBATCH --mail-user=jflores@unc.edu

set -euo pipefail

# =============================================================================
# Parameters
# =============================================================================

PROJECT_DIR="/work/users/j/p/jpflores/projects/STRS"
HIC_DIR="${PROJECT_DIR}/data/processed/hic/maps"
COOL_DIR="${PROJECT_DIR}/data/processed/hic/cool"
INSULATION_DIR="${PROJECT_DIR}/data/processed/hic/insulation"
EIGENVECTOR_DIR="${PROJECT_DIR}/data/processed/hic/eigenvectors"
FASTA="/proj/phanstiel_lab/Reference/human/hg38/fasta/GRCh38.primary_assembly.genome.fa"

CONTROL_HIC="${HIC_DIR}/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap_inter_30.hic"
SORBITOL_HIC="${HIC_DIR}/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap_inter_30.hic"

CONTROL_COOL="${COOL_DIR}/YAPP_HEK293_eGFP-YAP_Cai_control_megaMap.cool"
SORBITOL_COOL="${COOL_DIR}/YAPP_HEK293_eGFP-YAP_Cai_sorbitol_megaMap.cool"

RESOLUTION=10000
WINDOW=250000
N_SADDLE_BINS=50

# =============================================================================
# Environment
# =============================================================================

echo "[$(date)] Loading anaconda and activating cooltools_env..."
module load anaconda/2024.02
source activate cooltools_env

# =============================================================================
# Create output directories
# =============================================================================

mkdir -p "${COOL_DIR}" "${INSULATION_DIR}" "${EIGENVECTOR_DIR}"
mkdir -p "${PROJECT_DIR}/scripts/analysis/logs"

# =============================================================================
# Step 1: Convert .hic -> .cool
# =============================================================================
# hic2cool with -r produces a flat single-resolution .cool file.
# We delete any existing files first to avoid partial-write corruption.

echo "[$(date)] Converting control .hic to .cool..."
rm -f "${CONTROL_COOL}"
hic2cool convert "${CONTROL_HIC}" "${CONTROL_COOL}" -r "${RESOLUTION}"

echo "[$(date)] Converting sorbitol .hic to .cool..."
rm -f "${SORBITOL_COOL}"
hic2cool convert "${SORBITOL_HIC}" "${SORBITOL_COOL}" -r "${RESOLUTION}"

# =============================================================================
# Step 2: ICE balance both .cool files
# =============================================================================
# Required for cooltools insulation, eigs-cis, and saddle.
# --force overwrites any existing weight columns.

echo "[$(date)] Balancing control .cool..."
cooler balance --nproc "${SLURM_CPUS_PER_TASK}" --force "${CONTROL_COOL}"

echo "[$(date)] Balancing sorbitol .cool..."
cooler balance --nproc "${SLURM_CPUS_PER_TASK}" --force "${SORBITOL_COOL}"

# =============================================================================
# Step 3: Generate view file (autosomes only, in .cool chromosome order)
# =============================================================================
# The view file must match the exact chromosome order stored in the .cool file,
# otherwise cooltools raises a "not compatible viewframe" error.
# Order from cooler dump --table chroms:
#   chr1-7, chr8 comes after chrX in the file but we exclude chrX here,
#   then chr9, chr11, chr10, chr12-18, chr20, chr19, chrY (excluded), chr22, chr21
# We include autosomes only in exactly the order they appear in the .cool.

VIEW_FILE="${COOL_DIR}/hg38_autosomes.view.tsv"

echo "[$(date)] Generating autosome view file (matching .cool chromosome order)..."
printf \
"chr1\t0\t248956422\tchr1\n\
chr2\t0\t242193529\tchr2\n\
chr3\t0\t198295559\tchr3\n\
chr4\t0\t190214555\tchr4\n\
chr5\t0\t181538259\tchr5\n\
chr6\t0\t170805979\tchr6\n\
chr7\t0\t159345973\tchr7\n\
chr8\t0\t145138636\tchr8\n\
chr9\t0\t138394717\tchr9\n\
chr11\t0\t135086622\tchr11\n\
chr10\t0\t133797422\tchr10\n\
chr12\t0\t133275309\tchr12\n\
chr13\t0\t114364328\tchr13\n\
chr14\t0\t107043718\tchr14\n\
chr15\t0\t101991189\tchr15\n\
chr16\t0\t90338345\tchr16\n\
chr17\t0\t83257441\tchr17\n\
chr18\t0\t80373285\tchr18\n\
chr20\t0\t64444167\tchr20\n\
chr19\t0\t58617616\tchr19\n\
chr22\t0\t50818468\tchr22\n\
chr21\t0\t46709983\tchr21\n" > "${VIEW_FILE}"

# =============================================================================
# Step 4: Insulation score + boundary calls
# =============================================================================
# Diamond insulation score at 250kb window using ICE-balanced weights.
# --view restricts to autosomes in the correct order.
# --ignore-diags 2 skips self-ligation artifact bins near the diagonal.

echo "[$(date)] Computing insulation scores — control..."
cooltools insulation \
"${CONTROL_COOL}" "${WINDOW}" \
--view "${VIEW_FILE}" \
--ignore-diags 2 \
-p "${SLURM_CPUS_PER_TASK}" \
-o "${INSULATION_DIR}/control_insulation_${WINDOW}bp.tsv"

echo "[$(date)] Computing insulation scores — sorbitol..."
cooltools insulation \
"${SORBITOL_COOL}" "${WINDOW}" \
--view "${VIEW_FILE}" \
--ignore-diags 2 \
-p "${SLURM_CPUS_PER_TASK}" \
-o "${INSULATION_DIR}/sorbitol_insulation_${WINDOW}bp.tsv"

# =============================================================================
# Step 5: Generate bins BED and GC content track for PC1 orientation
# =============================================================================
# cooler dump extracts bin coordinates directly from the .cool file —
# guaranteed to match the exact bins used in eigs-cis.
# cooltools genome gc computes GC fraction per bin from the reference FASTA,
# used to orient PC1 so that positive = A compartment (GC-rich/active).

echo "[$(date)] Extracting genome bins from .cool (autosomes only)..."
echo -e "chrom\tstart\tend" \
> "${EIGENVECTOR_DIR}/genome_bins_${RESOLUTION}bp.bed"
cooler dump --table bins "${CONTROL_COOL}" \
| grep -P "^chr([1-9]|1[0-9]|2[012])\t" \
| cut -f1-3 \
>> "${EIGENVECTOR_DIR}/genome_bins_${RESOLUTION}bp.bed"

echo "[$(date)] Computing GC content per bin..."
cooltools genome gc \
"${EIGENVECTOR_DIR}/genome_bins_${RESOLUTION}bp.bed" \
"${FASTA}" \
> "${EIGENVECTOR_DIR}/gc_content_${RESOLUTION}bp.tsv"

# =============================================================================
# Step 6: Compute expected cis contact frequency (required for saddle)
# =============================================================================
# cooltools saddle requires a precomputed expected (average contact frequency
# by genomic distance) to normalize the observed contacts into O/E values.
# expected-cis computes this per chromosome and outputs a TSV.

echo "[$(date)] Computing expected cis — control..."
cooltools expected-cis \
"${CONTROL_COOL}" \
--view "${VIEW_FILE}" \
-p "${SLURM_CPUS_PER_TASK}" \
-o "${EIGENVECTOR_DIR}/control.expected.tsv"

echo "[$(date)] Computing expected cis — sorbitol..."
cooltools expected-cis \
"${SORBITOL_COOL}" \
--view "${VIEW_FILE}" \
-p "${SLURM_CPUS_PER_TASK}" \
-o "${EIGENVECTOR_DIR}/sorbitol.expected.tsv"

# =============================================================================
# Step 7: Compartment eigenvectors (PC1)
# =============================================================================
# eigs-cis performs PCA on the ICE-balanced observed/expected matrix per
# chromosome. --phasing-track orients PC1: positive = A (GC-rich/active).

echo "[$(date)] Computing eigenvectors — control..."
cooltools eigs-cis \
"${CONTROL_COOL}" \
--phasing-track "${EIGENVECTOR_DIR}/gc_content_${RESOLUTION}bp.tsv" \
--view "${VIEW_FILE}" \
--n-eigs 3 \
--out-prefix "${EIGENVECTOR_DIR}/control"

echo "[$(date)] Computing eigenvectors — sorbitol..."
cooltools eigs-cis \
"${SORBITOL_COOL}" \
--phasing-track "${EIGENVECTOR_DIR}/gc_content_${RESOLUTION}bp.tsv" \
--view "${VIEW_FILE}" \
--n-eigs 3 \
--out-prefix "${EIGENVECTOR_DIR}/sorbitol"

# =============================================================================
# Step 8: Saddle plots (compartment strength)
# =============================================================================
# cooltools saddle requires: the .cool, the eigenvector track, and the
# precomputed expected. We use control eigenvectors as reference for both
# conditions so quantile bins are directly comparable.

echo "[$(date)] Computing saddle plot — control..."
cooltools saddle \
"${CONTROL_COOL}" \
"${EIGENVECTOR_DIR}/control.cis.vecs.tsv::E1" \
"${EIGENVECTOR_DIR}/control.expected.tsv::balanced.avg" \
--view "${VIEW_FILE}" \
--n-bins "${N_SADDLE_BINS}" \
--qrange 0.02 0.98 \
--out-prefix "${EIGENVECTOR_DIR}/saddle_control"

echo "[$(date)] Computing saddle plot — sorbitol (control eigenvectors as reference)..."
cooltools saddle \
"${SORBITOL_COOL}" \
"${EIGENVECTOR_DIR}/control.cis.vecs.tsv::E1" \
"${EIGENVECTOR_DIR}/sorbitol.expected.tsv::balanced.avg" \
--view "${VIEW_FILE}" \
--n-bins "${N_SADDLE_BINS}" \
--qrange 0.02 0.98 \
--out-prefix "${EIGENVECTOR_DIR}/saddle_sorbitol"

# Extract saddledata matrix from .npz to TSV for import into R
echo "[$(date)] Extracting saddle matrices to TSV..."
python -c "
import numpy as np, os
d = '${EIGENVECTOR_DIR}'
for cond in ['control', 'sorbitol']:
    mat = np.load(os.path.join(d, f'saddle_{cond}.saddledump.npz'))['saddledata']
    np.savetxt(os.path.join(d, f'saddle_{cond}.saddledata.tsv'), mat, delimiter='\t')
    print(f'{cond}: shape {mat.shape}')
"

echo "[$(date)] All steps complete."
echo "Outputs:"
echo "  Insulation scores : ${INSULATION_DIR}/"
echo "  Expected cis      : ${EIGENVECTOR_DIR}/*.expected.tsv"
echo "  Eigenvectors      : ${EIGENVECTOR_DIR}/control.cis.vecs.tsv"
echo "                      ${EIGENVECTOR_DIR}/sorbitol.cis.vecs.tsv"
echo "  Saddle plots      : ${EIGENVECTOR_DIR}/saddle_control.*"
echo "                      ${EIGENVECTOR_DIR}/saddle_sorbitol.*"