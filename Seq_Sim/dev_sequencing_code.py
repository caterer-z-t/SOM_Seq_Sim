import pandas as pd
import numpy as np
import argparse
import os
from utils import DataPlotter, generate_synthetic_data, load_config


def main():
    # Argument parser to take either YAML config or individual arguments
    parser = argparse.ArgumentParser(
        description="Generate synthetic data and/or plot distributions."
    )

    # Add option for YAML config file
    parser.add_argument(
        "--config", help="Path to the configuration file (YAML format)."
    )

    # Add options for individual command-line arguments
    parser.add_argument("--input", help="Path to the original dataset (CSV file).")
    parser.add_argument(
        "--output", help="Path to store the generated synthetic data (CSV file)."
    )
    parser.add_argument(
        "--rows", type=int, help="Number of rows for the synthetic data."
    )
    parser.add_argument("--output-dir", help="Directory to save comparison plots.")

    args = parser.parse_args()

    # Ensure that either config or input argument is provided
    if not args.config and not args.input:
        parser.error("You must provide either --config or --input.")

    # Load configuration from YAML or individual arguments
    input_path, output_path, rows, output_dir = load_config(args)

    if rows and not output_path:
        raise ValueError("The 'output' path must be specified if 'rows' is specified.")

    # Read the input file
    df = pd.read_csv(input_path, sep="\t")

    synthetic_df = None
    if rows:
        # Generate synthetic data
        synthetic_df = generate_synthetic_data(df, rows)

        # Save the synthetic data to output file
        if output_path:
            # Ensure output directory exists
            output_dir_path = os.path.dirname(output_path)
            if not os.path.exists(output_dir_path) and output_dir_path:
                os.makedirs(output_dir_path)
            synthetic_df.to_csv(output_path, index=False)

    # Create an instance of DataPlotter
    plotter = DataPlotter(metadata=df, synthetic_metadata=synthetic_df)

    # Plot numerical distributions
    plotter.plot_numerical_distributors()

    # Get a random row from the dataframe
    random_row = df.iloc[np.random.randint(0, df.shape[0])]

    plotter_with_rows = DataPlotter(metadata=df, target_row=random_row)

    # Plot numerical distributions with a target row
    plotter_with_rows.plot_numerical_distributors()

    # Plot comparisons if synthetic data was generated
    if synthetic_df is not None and output_dir:
        plotter.plot_comparisons(output_dir=output_dir)


if __name__ == "__main__":
    main()
