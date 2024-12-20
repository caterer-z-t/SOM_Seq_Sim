# **Sequence Simulation**

The `Seq_Sim` directory contains the code and configuration necessary for generating and saving synthetic biological datasets. This process involves simulating different conditions such as sample size and fold change, based on parameters defined in a YAML configuration file. The Python script `seq_sim.py` orchestrates the data generation, while a shell script `simulation.sh` automates the execution of the simulation using parameters from the configuration file.

## **Table of Contents**
- [Directory Structure](#directory-structure)
- [File Descriptions](#file-descriptions)
- [How it Works](#how-it-works)
- [How to Use](#how-to-use)
- [Customization](#customizing-the-simulation)
- [Conclusions](#conclusion)


## Directory Structure

Here’s an overview of the key files and their roles:

``` bash
Seq_Sim/
├── utils/
│   └── seq_sim_utils.py        # Utility functions for handling configuration, validation, and feature generation
├── zhang_lab/                  # Code written by the Zhang Lab, the python files in this parent directory were modified from these scripts
│   ├── function.R              # Utitlity functions for generating the single cell sequencing
│   ├── simulation.R            # Main R script to call the utility function and create the sequencing datasets and figures
│   ├── simulation.sh           # Shell script to run simulations with various parameters
│   └── config.yml              # Configuration file containing simulation parameters
├── simulation.sh               # Shell script to run simulations with various parameters
├── config.yml                  # Configuration file containing simulation parameters
└── seq_sim.py                  # Main Python script for data simulation
```

### File Descriptions

1. **`seq_sim.py`**:
   - This Python script is the main entry point for the simulation. It handles loading configuration files, validating arguments, and calling functions to generate the synthetic dataset based on user-specified parameters (e.g., number of samples, fold change).
     - **Main functions used from `utils/seq_sim_utils.py`**:
     - `arg_parser()`: Parses command-line arguments.
     - `validate_arguments()`: Validates the provided arguments.
     - `load_config()`: Loads the YAML configuration file.
     - `generate_and_save_features()`: Simulates data based on the parameters and saves the output.

2. **`utils/seq_sim_utils.py`**:
   - This utility script contains helper functions for reading the configuration file, validating inputs, and generating data. The functions make it easier to modularize the code and keep the main script clean and readable.

3. **`simulation.sh`**:
   - This is a shell script that automates the running of the simulation. It reads parameters from the YAML configuration file and executes the Python simulation script with different combinations of `num_samples`, `fold_changes`, and `config.yml`.
   - The script also handles logging by writing the simulation results to a log file and ensures that the necessary file permissions are set.

4. **`config.yml`**:
   - This YAML configuration file defines the parameters for the simulation. It specifies the number of samples, fold changes, paths to data and scripts, and other simulation-specific settings like log file paths and dataset parameters.
   - **Important sections**:
     - `num_samples`: A list of integers specifying how many samples to simulate.
     - `fold_changes`: A list of fold change values that will be used in the simulation.
     - `log_file`: Path to the log file where errors and output are saved.
     - `file_path_to_simulation`: Path to the `seq_sim.py` script, relative to the `Seq_Sim/` directory.

   - Example of the structure and layout of this configuration file is shown below. 
   ``` yaml
   log_file: "error.log"
   data_file_path: "./data/"
   file_path_to_simulation: "Seq_Sim/seq_sim.py"
   functions_script_path: "Seq_Sim/utils/seq_sim_utils.py"
   file_prefix: "sim_data"

   num_samples:
   - 10
   - 20
   - 30
   fold_changes:
   - 0.1
   - 0.75
   - 1.5
   - 3

   n_threads: 4

   dummy_dataset_params:
   n_cells: 100 # cells per individual #PER CELLTYPE
   sd_celltypes: 0.1 # standard deviation of number of cells
   n_major_cell_types: 7
   n_minor_cell_types: 3
   relative_abundance: 0.4 # ratio between major and rare
   n_major_diff_celltypes: 1
   n_minor_diff_celltypes: 1
   n_batchs: 4
   prop_sex: 0.5 # proportion of feature 1 if feature 1 is categorical variable
   prop_disease: 0.5 # proportion of feature 2 if feature 2 is categorical variable
   seed: 1234
   n_features: 1000

   variance_attributes: 
   cluster_ratio: 0.7

   ratio_variance: 0.1

   column_information:
   cluster_col: "cell_type"
   disease_col: "disease"
   individual_col: "subject_id"

   files_to_save:
   feature_matrix: True
   latent_factors: True
   ```

## How It Works

1. **Running the Simulation**:
   The simulation is executed by running the `simulation.sh` script. This script:
   - Reads the `config.yml` file to get the list of `num_samples`, `fold_changes`, and other parameters.
   - For each combination of `num_samples` and `fold_change`, it runs the `seq_sim.py` Python script, passing the parameters as command-line arguments.
   - The simulation results are saved in the specified log file and other output files (e.g., feature matrices, latent factors).

2. **Data Generation**:
   In `seq_sim.py`, the main function orchestrates the data generation. It first validates the arguments, loads the configuration, and calls the `generate_and_save_features()` function to simulate synthetic data based on the loaded configuration and parameters.

3. **Logging and Output**:
   The `simulation.sh` script ensures that the log file is created (if it doesn't already exist) and that it has the correct permissions. All simulation output, including error messages and progress, is logged in the specified log file.

## How to Use

### Prerequisites

1. **Install Dependencies**:
   You need Python 3.12 (or compatible version) and some dependencies installed via Conda. You can create the Conda environment and install dependencies by running:

   ```bash
   conda env create -f environment.yml
   conda activate som-sim-env
   ```

2. **Configuration**:
   The `config.yml` file contains all the necessary parameters for the simulation, including the path to the simulation Python script, the number of samples, fold changes, and more. You can modify this file to adjust the simulation settings. See [File Descriptions](#file-descriptions) or [config.yml](config.yml) for more information. 

3. **Run the Simulation**:
   After ensuring that the Conda environment is activated and the configuration is set, you can run the simulation by executing the `simulation.sh` script:

   ```bash
   ./simulation.sh --config <path_to_config_file>
   ```

   This will run the simulation using the parameters specified in the `config.yml` file.

   Alternatively, you can run a single simulation with a single value for `num_samples` and `fold_change` with:
   ``` bash 
   python seq_sim.py \ 
        --num_samples 30                   \ # num samples 
        --fold_change 0.5                  \ # fold change between disease and healthy samples 
        --config_file <path_to_config_file>    # configuration file
    ```

### Customizing the Simulation

- **Changing Simulation Parameters**:
  To change the number of samples or fold changes, simply modify the `num_samples` and `fold_changes` arrays in the `config.yml` file. These arrays define the values that will be used for each simulation run.

- **Adding/Modifying Dataset Parameters**:
  If you need to adjust how the synthetic dataset is generated (e.g., number of cell types, relative abundance), modify the `dummy_dataset_params` section in `config.yml`.

## Conclusion

This directory provides a streamlined way to simulate synthetic biological datasets with different parameters. By modifying the configuration file, you can easily adjust the simulation settings to meet your needs. The use of Conda ensures that all dependencies are managed, and the shell script automates the process of running simulations with varying parameters.
