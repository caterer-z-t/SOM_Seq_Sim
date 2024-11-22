#!/bin/bash

# Function to display help message
show_help() {
    echo "Usage: ./your_script.sh [options]"
    echo
    echo "Options:"
    echo "  -c, --config <path>   Specify the path to the YAML configuration file."
    echo "  -h, --help            Display this help message."
}

# Default values
config_file="config.yml"

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -c|--config) config_file="$2"; shift ;;
        -h|--help) show_help; exit 0 ;;
        *) echo "Unknown parameter: $1"; show_help; exit 1 ;;
    esac
    shift
done

# Exit immediately if a command exits with a non-zero status
set -e

# Load parameters from the YAML configuration file
eval "$(yq e '.environment' $config_file) env"
eval "$(yq e '.log_file' $config_file) log_file"
eval "$(yq e '.data_file_path' $config_file) data_file_path"
eval "$(yq e '.file_path_to_simulation' $r_script_file) r_script_file"
num_samples=($(yq e '.num_samples[]' $config_file))
fold_changes=($(yq e '.fold_changes[]' $config_file))

# Activate the Conda environment
eval "$(conda shell.bash hook)"
conda activate "$environment"

# Redirect output and errors to a log file
exec &> "$log_file"

# Loop through the arrays and generate datasets
for num_samp in "${num_samples[@]}"; do
    echo "Number of samples: $num_samp"
    for fold_change in "${fold_changes[@]}"; do
        echo "Fold change: $fold_change"
        
        # Construct the filenames based on the parameters
        file_prefix="simulated_data_numSamp${num_samp}_foldChange${fold_change}"
        
        # Run the R script with additional arguments for filenames
        Rscript "$r_script_file" "$num_samp" "$fold_change" "$data_file_path" "$file_prefix"
    done
done
