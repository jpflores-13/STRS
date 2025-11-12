#!/bin/bash
#SBATCH --job-name=homer_top500
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=homer_top500_%j.out
#SBATCH --error=homer_top500_%j.err

# Load HOMER
module load homer/5.1

# Set working directory
cd /work/users/j/p/jpflores/projects/STRS

# Create output directories
mkdir -p data/processed/cutntag/homer_motifs/top500
mkdir -p data/processed/cutntag/homer_preparsed

echo "========================================="
echo "HOMER Motif Analysis: Top 500 Peaks"
echo "========================================="
echo ""

# Array of samples to process
proteins=(ctcf h3k27ac rad21 yap1)
conditions=(control sorbitol)

# Loop through each sample
for protein in "${proteins[@]}"; do
  for cond in "${conditions[@]}"; do
    
    sample="${protein}_${cond}"
    bed_file="data/processed/cutntag/homer_input/top500_${sample}.bed"
    output_dir="data/processed/cutntag/homer_motifs/top500/${sample}"
    
    # Check if BED file exists
    if [ ! -f "$bed_file" ]; then
      echo "WARNING: File not found: $bed_file"
      continue
    fi
    
    echo "========================================="
    echo "Processing: ${sample}"
    echo "========================================="
    
    # Create output directory
    mkdir -p "$output_dir"
    
    # Run HOMER
    findMotifsGenome.pl \
      "$bed_file" \
      hg38 \
      "$output_dir" \
      -size 200 \
      -mask \
      -p 8 \
      -preparsedDir data/processed/cutntag/homer_preparsed
    
    echo "${sample} complete!"
    echo ""
    
  done
done

echo "========================================="
echo "Summary"
echo "========================================="
echo "All results saved to: data/processed/cutntag/homer_motifs/top500/"
echo ""
echo "To view results, open the homerResults.html files in each directory"
echo ""
echo "Job complete!"
