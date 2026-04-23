#!/bin/bash

#SBATCH --job-name=domain_analysis_YAPdTAD
#SBATCH --output=/work/users/j/p/jpflores/projects/STRS/scripts/analysis/logs/domain_analysis_YAPdTAD_%j.out
#SBATCH --error=/work/users/j/p/jpflores/projects/STRS/scripts/analysis/logs/domain_analysis_YAPdTAD_%j.err
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
INSULATION_DIR="${PROJECT_DIR}/data/processed/hic/insulation_YAPdTAD"
EIGENVECTOR_DIR="${PROJECT_DIR}/data/processed/hic/eigenvectors_YAPdTAD"
FASTA="/proj/phanstiel_lab/Reference/human/hg38/fasta/GRCh38.primary_assembly.genome.fa"

# Single-rep YAPdTAD .hic files (no megaMap available)
CONTROL_HIC="${HIC_DIR}/YAPP_HEK293_eGFP-YAPdTAD_Cai_control_1_1_inter_30.hic"
SORBITOL_HIC="${HIC_DIR}/YAPP_HEK293_eGFP-YAPdTAD_Cai_sorbitol_1_1_inter_30.hic"

CONTROL_COOL="${COOL_DIR}/YAPP_HEK293_eGFP-YAPdTAD_Cai_control.cool"
SORBITOL_COOL="${COOL_DIR}/YAPP_HEK293_eGFP-YAPdTAD_Cai_sorbitol.cool"

# Reuse the view file and GC content track already generated for eGFP-YAP
# (same genome, same resolution — no need to regenerate)
VIEW_FILE="${COOL_DIR}/hg38_autosomes.view.tsv"
GC_TRACK="${PROJECT_DIR}/data/processed/hic/eigenvectors/gc_content_10000bp.tsv"

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

echo "[$(date)] Converting YAPdTAD control .hic to .cool..."
rm -f "${CONTROL_COOL}"
hic2cool convert "${CONTROL_HIC}" "${CONTROL_COOL}" -r "${RESOLUTION}"

echo "[$(date)] Converting YAPdTAD sorbitol .hic to .cool..."
rm -f "${SORBITOL_COOL}"
hic2cool convert "${SORBITOL_HIC}" "${SORBITOL_COOL}" -r "${RESOLUTION}"

# =============================================================================
# Step 2: ICE balance both .cool files
# =============================================================================

echo "[$(date)] Balancing YAPdTAD control .cool..."
cooler balance --nproc "${SLURM_CPUS_PER_TASK}" --force "${CONTROL_COOL}"

echo "[$(date)] Balancing YAPdTAD sorbitol .cool..."
cooler balance --nproc "${SLURM_CPUS_PER_TASK}" --force "${SORBITOL_COOL}"

# =============================================================================
# Step 3: Insulation scores + boundary calls
# =============================================================================

echo "[$(date)] Computing insulation scores -- YAPdTAD control..."
cooltools insulation \
  "${CONTROL_COOL}" "${WINDOW}" \
  --view "${VIEW_FILE}" \
  --ignore-diags 2 \
  -p "${SLURM_CPUS_PER_TASK}" \
  -o "${INSULATION_DIR}/control_insulation_${WINDOW}bp.tsv"

echo "[$(date)] Computing insulation scores -- YAPdTAD sorbitol..."
cooltools insulation \
  "${SORBITOL_COOL}" "${WINDOW}" \
  --view "${VIEW_FILE}" \
  --ignore-diags 2 \
  -p "${SLURM_CPUS_PER_TASK}" \
  -o "${INSULATION_DIR}/sorbitol_insulation_${WINDOW}bp.tsv"

# =============================================================================
# Step 4: Compute expected cis contact frequency
# =============================================================================

echo "[$(date)] Computing expected cis -- YAPdTAD control..."
cooltools expected-cis \
  "${CONTROL_COOL}" \
  --view "${VIEW_FILE}" \
  -p "${SLURM_CPUS_PER_TASK}" \
  -o "${EIGENVECTOR_DIR}/control.expected.tsv"

echo "[$(date)] Computing expected cis -- YAPdTAD sorbitol..."
cooltools expected-cis \
  "${SORBITOL_COOL}" \
  --view "${VIEW_FILE}" \
  -p "${SLURM_CPUS_PER_TASK}" \
  -o "${EIGENVECTOR_DIR}/sorbitol.expected.tsv"

# =============================================================================
# Step 5: Compartment eigenvectors (PC1)
# =============================================================================

echo "[$(date)] Computing eigenvectors -- YAPdTAD control..."
cooltools eigs-cis \
  "${CONTROL_COOL}" \
  --phasing-track "${GC_TRACK}" \
  --view "${VIEW_FILE}" \
  --n-eigs 3 \
  --out-prefix "${EIGENVECTOR_DIR}/control"

echo "[$(date)] Computing eigenvectors -- YAPdTAD sorbitol..."
cooltools eigs-cis \
  "${SORBITOL_COOL}" \
  --phasing-track "${GC_TRACK}" \
  --view "${VIEW_FILE}" \
  --n-eigs 3 \
  --out-prefix "${EIGENVECTOR_DIR}/sorbitol"

# =============================================================================
# Step 6: Saddle plots (compartment strength)
# =============================================================================

echo "[$(date)] Computing saddle -- YAPdTAD control..."
cooltools saddle \
  "${CONTROL_COOL}" \
  "${EIGENVECTOR_DIR}/control.cis.vecs.tsv::E1" \
  "${EIGENVECTOR_DIR}/control.expected.tsv::balanced.avg" \
  --view "${VIEW_FILE}" \
  --n-bins "${N_SADDLE_BINS}" \
  --qrange 0.02 0.98 \
  --out-prefix "${EIGENVECTOR_DIR}/saddle_control"

echo "[$(date)] Computing saddle -- YAPdTAD sorbitol (control eigenvectors as reference)..."
cooltools saddle \
  "${SORBITOL_COOL}" \
  "${EIGENVECTOR_DIR}/control.cis.vecs.tsv::E1" \
  "${EIGENVECTOR_DIR}/sorbitol.expected.tsv::balanced.avg" \
  --view "${VIEW_FILE}" \
  --n-bins "${N_SADDLE_BINS}" \
  --qrange 0.02 0.98 \
  --out-prefix "${EIGENVECTOR_DIR}/saddle_sorbitol"

# =============================================================================
# Step 7: Extract saddle matrices from .npz to TSV for R
# =============================================================================

echo "[$(date)] Extracting saddle matrices to TSV..."
for COND in control sorbitol; do
    python3 - "${EIGENVECTOR_DIR}" "${COND}" << 'PYEOF'
import sys, numpy as np, os
d, cond = sys.argv[1], sys.argv[2]
npz_path = os.path.join(d, f"saddle_{cond}.saddledump.npz")
mat = np.load(npz_path)["saddledata"]
out_path = os.path.join(d, f"saddle_{cond}.saddledata.tsv")
np.savetxt(out_path, mat, delimiter="\t")
print(f"{cond}: shape {mat.shape}, saved to {out_path}")
PYEOF
done

echo "[$(date)] All steps complete."
echo "Outputs:"
echo "  Insulation scores : ${INSULATION_DIR}/"
echo "  Expected cis      : ${EIGENVECTOR_DIR}/*.expected.tsv"
echo "  Eigenvectors      : ${EIGENVECTOR_DIR}/control.cis.vecs.tsv"
echo "                      ${EIGENVECTOR_DIR}/sorbitol.cis.vecs.tsv"
echo "  Saddle data       : ${EIGENVECTOR_DIR}/saddle_control.saddledata.tsv"
echo "                      ${EIGENVECTOR_DIR}/saddle_sorbitol.saddledata.tsv"