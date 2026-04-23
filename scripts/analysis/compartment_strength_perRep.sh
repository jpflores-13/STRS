#!/bin/bash

#SBATCH --job-name=compartment_strength_perRep
#SBATCH --output=/work/users/j/p/jpflores/projects/STRS/scripts/analysis/logs/compartmentStrength_perRep_%j.out
#SBATCH --error=/work/users/j/p/jpflores/projects/STRS/scripts/analysis/logs/compartmentStrength_perRep_%j.err
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
OUT_DIR="${PROJECT_DIR}/data/processed/hic/compartment_strength"

## Reuse view file, GC track, and eigenvectors from the megaMap eGFP-YAP run
VIEW_FILE="${COOL_DIR}/hg38_autosomes.view.tsv"
GC_TRACK="${PROJECT_DIR}/data/processed/hic/eigenvectors/gc_content_10000bp.tsv"
## Control eigenvectors used as shared reference for all saddle calls so
## quantile bins are directly comparable across samples
EIG_REF="${PROJECT_DIR}/data/processed/hic/eigenvectors/control.cis.vecs.tsv"

RESOLUTION=10000
N_SADDLE_BINS=50

# =============================================================================
# Environment
# =============================================================================

echo "[$(date)] Loading anaconda and activating cooltools_env..."
module load anaconda/2024.02
source activate cooltools_env

mkdir -p "${COOL_DIR}" "${OUT_DIR}"
mkdir -p "${PROJECT_DIR}/scripts/analysis/logs"

# =============================================================================
# Define per-replicate samples
# Format: "LABEL:HIC_FILE"
# =============================================================================

declare -a SAMPLES=(
  "eGFP-YAP_ctrl_rep1:YAPP_HEK293_eGFP-YAP_Cai_control_1_2_inter_30.hic"
  "eGFP-YAP_ctrl_rep2:YAPP_HEK293_eGFP-YAP_Cai_control_2_2_inter_30.hic"
  "eGFP-YAP_ctrl_rep3:YAPP_HEK293_eGFP-YAP_Cai_control_3_2_inter_30.hic"
  "eGFP-YAP_ctrl_rep4:YAPP_HEK293_eGFP-YAP_Cai_control_4_1_inter_30.hic"
  "eGFP-YAP_sorb_rep1:YAPP_HEK293_eGFP-YAP_Cai_sorbitol_5_2_inter_30.hic"
  "eGFP-YAP_sorb_rep2:YAPP_HEK293_eGFP-YAP_Cai_sorbitol_6_2_inter_30.hic"
  "eGFP-YAP_sorb_rep3:YAPP_HEK293_eGFP-YAP_Cai_sorbitol_7_2_inter_30.hic"
  "eGFP-YAP_sorb_rep4:YAPP_HEK293_eGFP-YAP_Cai_sorbitol_8_1_inter_30.hic"
  "YAPdTAD_ctrl_rep1:YAPP_HEK293_eGFP-YAPdTAD_Cai_control_1_1_inter_30.hic"
  "YAPdTAD_sorb_rep1:YAPP_HEK293_eGFP-YAPdTAD_Cai_sorbitol_1_1_inter_30.hic"
)

# =============================================================================
# Process each replicate
# =============================================================================

for ENTRY in "${SAMPLES[@]}"; do

  LABEL="${ENTRY%%:*}"
  HIC_FILE="${HIC_DIR}/${ENTRY##*:}"
  COOL_FILE="${COOL_DIR}/${LABEL}.cool"
  EXPECTED_FILE="${OUT_DIR}/${LABEL}.expected.tsv"
  SADDLE_PREFIX="${OUT_DIR}/saddle_${LABEL}"
  SADDLE_TSV="${OUT_DIR}/saddle_${LABEL}.saddledata.tsv"

  echo "[$(date)] Processing: ${LABEL}"

  # Step 1: Convert .hic -> .cool
  echo "  Converting to .cool..."
  rm -f "${COOL_FILE}"
  hic2cool convert "${HIC_FILE}" "${COOL_FILE}" -r "${RESOLUTION}"

  # Step 2: ICE balance
  echo "  Balancing..."
  cooler balance --nproc "${SLURM_CPUS_PER_TASK}" --force "${COOL_FILE}"

  # Step 3: Expected cis
  echo "  Computing expected cis..."
  cooltools expected-cis \
    "${COOL_FILE}" \
    --view "${VIEW_FILE}" \
    -p "${SLURM_CPUS_PER_TASK}" \
    -o "${EXPECTED_FILE}"

  # Step 4: Saddle using shared control eigenvectors as reference
  echo "  Computing saddle..."
  cooltools saddle \
    "${COOL_FILE}" \
    "${EIG_REF}::E1" \
    "${EXPECTED_FILE}::balanced.avg" \
    --view "${VIEW_FILE}" \
    --n-bins "${N_SADDLE_BINS}" \
    --qrange 0.02 0.98 \
    --out-prefix "${SADDLE_PREFIX}"

  # Step 5: Extract saddledata matrix to TSV for R
  echo "  Extracting saddle matrix to TSV..."
  python3 - "${SADDLE_PREFIX}.saddledump.npz" "${SADDLE_TSV}" << 'PYEOF'
import sys, numpy as np
npz_path, out_path = sys.argv[1], sys.argv[2]
mat = np.load(npz_path)["saddledata"]
np.savetxt(out_path, mat, delimiter="\t")
print(f"Shape {mat.shape} -> {out_path}")
PYEOF

  echo "  Done: ${LABEL}"

done

echo "[$(date)] All replicates complete."
echo "Saddle TSVs in: ${OUT_DIR}/"