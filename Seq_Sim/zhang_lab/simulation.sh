#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Declare arrays
declare -a num_samples=(10 20 30)
declare -a fold_changes=(0.1 0.75 1.5 3)

# Activate the Conda environment
eval "$(conda shell.bash hook)"
conda activate /opt/anaconda3/envs/zhang_lab

# Redirect output and errors to a log file
exec &> /Users/zc/Library/CloudStorage/OneDrive-UCB-O365/Zhang_lab/zac/Zac_SHAP/miloR/error_log.txt

# Loop through the arrays and generate datasets
for num_samp in "${num_samples[@]}"; do
    echo "Number of samples: $num_samp"
    for fold_change in "${fold_changes[@]}"; do
        echo "Fold change: $fold_change"
        # Run the R script
        Rscript /Users/zc/Library/CloudStorage/OneDrive-UCB-O365/Zhang_lab/jade/xcell_results/simulation/simulation_bench.R "$num_samp" "$fold_change"
    done
done

# Run the CellPhenoX Python script
# python3 ./xcell_functions/main.py
