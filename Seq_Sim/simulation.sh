#!/bin/bash

# Function to display help message
show_help() {
    echo "Usage: ./simulation.sh [options]"
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
log_file=$(yq e '.log_file' $config_file)
num_samples=($(yq e '.num_samples[]' $config_file))
fold_changes=($(yq e '.fold_changes[]' $config_file))
seq_sim=$(yq e '.file_path_to_simulation' $config_file)


# Ensure log file exists
if [[ ! -f "$log_file" ]]; then
    touch "$log_file"
fi

# Ensure log file has the correct permissions
chmod u+rw "$log_file"

# Redirect output and errors to a log file
exec &> "$log_file"

# Loop through the arrays and generate datasets
for num_samp in "${num_samples[@]}"; do
    echo "Number of samples: $num_samp"
    for fold_change in "${fold_changes[@]}"; do
        echo "Fold change: $fold_change"
        
        # Construct the filenames based on the parameters
        file_prefix="sim${num_samp}_fc${fold_change}"
        
        # Run the R script with additional arguments for filenames
        python "$seq_sim" "$num_samp" "$fold_change" "$config_file"
    done
done
