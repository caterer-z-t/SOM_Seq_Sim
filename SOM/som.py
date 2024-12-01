"""
SOM Hyperparameter Tuning and Fitting
=====================================

Description:
------------
This script provides a command-line interface for fitting and visualizing a Self-Organizing Map
(SOM). It supports both direct SOM fitting and hyperparameter tuning, depending on the inputs 
provided.

Inputs:
=======
- **Training Data (-t / --train_dat)** *(required)*:  
        Path to a CSV file containing the numerical training data. Each row represents a data
        point, and each column represents a feature.

- **Categorical Data (-c / --other_dat)** *(optional)*:  
        Path to a CSV file containing additional categorical or other metadata associated with
        the training data. If only one column is provided, it will still be treated as a DataFrame.

- **Output Directory (-o / --output_directory)** *(required)*:  
        Path to the directory where results, including metrics and visualizations, will be saved.

- **Scaling Method (-s / --scale_method)** *(required)*:  
        The method used to scale the data before training the SOM. Accepts one or more values:  
        `zscore` or `minmax`.

- **x-dimension (-x / --x_dim)** *(required)*:  
        Number of neurons in the x-dimension of the SOM grid. Accepts one or more integer values.

- **y-dimension (-y / --y_dim)** *(required)*:  
        Number of neurons in the y-dimension of the SOM grid. Accepts one or more integer values.

- **Topology (-p / --topology)** *(required)*:  
        The grid topology for the SOM. Accepts one or more values: `rectangular` or `hexagonal`.

- **Neighborhood Function (-n / --neighborhood_fnc)** *(required)*:  
        The neighborhood function for training the SOM. Accepts one or more values: `gaussian` or
        `bubble`.

- **Epochs (-e / --epochs)** *(required)*:  
        Number of training epochs for the SOM. Accepts one or more integer values.

- **Component Plane Plots (-m / --plot_component_planes)** *(optional)*:  
        Whether to generate component plane plots for each feature in the training data. Defaults
        to `True`.

Usage:
======
Run the script from the command line with the required input paths and hyperparameters. If a single
value is passed for each hyperparameter, the SOM is fitted directly. If multiple values are 
provided, the script performs hyperparameter tuning to find the optimal configuration.

Example command to perform hyperparameter tuning:
-------------------------------------------------
python som_script.py -t train_data.csv -c categorical_data.csv -o output_dir \
    -s zscore minmax -x 3 5 7 -y 2 4 -p rectangular hexagonal \
    -n gaussian bubble -e 50 100

Example command to fit a SOM with predetermined hyperparameters:
----------------------------------------------------------------
python som_script.py -t train_data.csv -c categorical_data.csv -o output_dir \
    -s zscore -x 5  -y 4 -p hexagonal -n gaussian bubble -e 50 100
"""

# Standard imports
import argparse
from itertools import product
import os
import sys

# Third party imports
import pandas as pd

# Local imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from SOM.utils.som_utils import SOM


def main():

    """
    Main function to fit and/or tune a Self-Organizing Map (SOM) using user-provided inputs.

    This function reads in training data and optional categorical data, accepts user-specified 
    hyperparameters, and performs SOM fitting or hyperparameter tuning based on the number of 
    hyperparameter values provided. Results, including performance metrics and visualizations,
    are saved to the specified output directory.  

    Functionality:
    --------------
    - If **single values** are provided for all hyperparameters:
        - The SOM is fitted directly with the specified parameters, and results are saved and 
          visualized.
    - If **multiple values** are provided for any hyperparameter:
        - Performs hyperparameter tuning using all combinations of provided values.
        - Evaluates each combination based on Percent Variance Explained (PVE) and Topographic 
          Error.
        - Identifies and outputs the best hyperparameter combination and retrains the SOM using 
          these parameters.

    Outputs:
    --------
    - **Metrics**:  
      Displays PVE and Topographic Error for the fitted SOM or for each hyperparameter combination
      during tuning.

    - **Best Parameters** (if tuning):  
      Outputs the combination of hyperparameters that yielded the best performance.

    - **Visualizations**:  
      - Component plane plots for each feature (if enabled).  
      - SOM grid mapping of categorical data (if provided).

    Notes:
    ------
    - Hyperparameter tuning may take considerable time for large datasets or extensive
      hyperparameter grids.
    """

    parser = argparse.ArgumentParser(
        description=(
            "Fit a Self-Organizing Map (SOM) with specified hyperparameters. Performs"
            "hyperparameter tuning if multiple values/options for hyperparameters are provided."
        )
    )

    # Required input paths
    parser.add_argument(
        '-t', '--train_dat',
        type=str,
        help='Path to the CSV file containing training data',
        required=True,
        metavar=""
    )
    parser.add_argument(
        '-c', '--other_dat',
        type=str,
        help='Path to the CSV file containing other/categorical data',
        required=False,
        metavar=""
    )
    parser.add_argument(
        '-o', '--output_directory',
        type=str,
        help='Path to the output directory for saving results',
        required=True,
        metavar=""
    )

    # Required hyperparameters
    parser.add_argument(
        '-s', '--scale_method',
        type=str,
        nargs='+',
        required=True,
        help="Scaling method for the data ('zscore' or 'minmax'). Specify one or more values."
    )
    parser.add_argument(
        '-x', '--x_dim',
        type=int,
        nargs='+',
        required=True,
        help="Number of neurons in the x-dimension. Specify one or more values."
    )
    parser.add_argument(
        '-y', '--y_dim',
        type=int,
        nargs='+',
        required=True,
        help="Number of neurons in the y-dimension. Specify one or more values."
    )
    parser.add_argument(
        '-p', '--topology',
        type=str,
        nargs='+',
        required=True,
        help="Grid topology ('rectangular' or 'hexagonal'). Specify one or more values."
    )
    parser.add_argument(
        '-n', '--neighborhood_fnc',
        type=str,
        nargs='+',
        required=True,
        help="Neighborhood function ('gaussian' or 'bubble'). Specify one or more values."
    )
    parser.add_argument(
        '-e', '--epochs',
        type=int,
        nargs='+',
        required=True,
        help="Number of training epochs. Specify one or more values."
    )
    parser.add_argument(
        '-m', '--plot_component_planes',
        type=str,
        default=True,
        help="Whether to create component plane plot for every feature in training data."
    )
    args = parser.parse_args()

    # Load the data
    train_dat = pd.read_csv(args.train_dat)

    if args.other_dat:
        other_dat = pd.read_csv(args.other_dat)
        # Ensure DataFrame structure even if only 1 column
        if isinstance(other_dat, pd.Series):
            other_dat = other_dat.to_frame()

    # Check if hyperparameter tuning is needed
    hyperparameter_combinations = list(product(
        args.scale_method,
        args.x_dim,
        args.y_dim,
        args.topology,
        args.neighborhood_fnc,
        args.epochs
    ))

    if len(hyperparameter_combinations) > 1:
        # Hyperparameter tuning
        best_params = None
        best_score = -float("inf")

        for params in hyperparameter_combinations:
            scale_method, x_dim, y_dim, topology, neighborhood_fnc, epochs = params

            # Train SOM with current parameters
            som = SOM(
                train_dat=train_dat,
                other_dat=other_dat,
                scale_method=scale_method,
                x_dim=x_dim,
                y_dim=y_dim,
                topology=topology,
                neighborhood_fnc=neighborhood_fnc,
                epochs=epochs
            )
            som.train_map()

            # Evaluate performance
            pve = som.calculate_percent_variance_explained()
            topographic_error = som.calculate_topographic_error()
            score = pve - topographic_error * 100

            print(f"Tested params: {params} | Score: {score:.2f}")

            # Update best parameters if score improves
            if score > best_score:
                best_score = score
                best_params = params

        print("\nBest Hyperparameters:")
        print(f"\tScale Method: {best_params[0]}")
        print(f"\tx_dim: {best_params[1]}")
        print(f"\ty_dim: {best_params[2]}")
        print(f"\tTopology: {best_params[3]}")
        print(f"\tNeighborhood Function: {best_params[4]}")
        print(f"\tEpochs: {best_params[5]}")
        print(f"\tBest Score: {best_score:.2f}")

        # Train SOM with best parameters
        best_som = SOM(
            train_dat=train_dat,
            other_dat=other_dat,
            scale_method=best_params[0],
            x_dim=best_params[1],
            y_dim=best_params[2],
            topology=best_params[3],
            neighborhood_fnc=best_params[4],
            epochs=best_params[5]
        )
        best_som.train_map()

        # Save and visualize results
        if args.plot_component_planes:
            best_som.plot_component_planes(
                output_dir=args.output_directory
            )
        if args.other_dat:
            best_som.plot_categorical_data(
                output_dir=args.output_directory
            )

    else:
        # Single-value parameters: Direct SOM fitting
        som = SOM(
            train_dat=train_dat,
            other_dat=other_dat,
            scale_method=args.scale_method[0],
            x_dim=args.x_dim[0],
            y_dim=args.y_dim[0],
            topology=args.topology[0],
            neighborhood_fnc=args.neighborhood_fnc[0],
            epochs=args.epochs[0]
        )
        som.train_map()

        # Save and visualize results
        if args.plot_component_planes:
            som.plot_component_planes(
                output_dir=args.output_directory
            )
        if args.other_dat:
            som.plot_categorical_data(
                output_dir=args.output_directory
            )

        # Output metrics
        pve = som.calculate_percent_variance_explained()
        topographic_error = som.calculate_topographic_error()

        print("\nSOM Performance:")
        print(f"Percent variance explained = {pve}%")
        print(f"Topographic error = {topographic_error * 100}%")


if __name__ == "__main__":
    main()
