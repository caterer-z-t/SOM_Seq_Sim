import pandas as pd
import numpy as np
import argparse
import os
from utils import DataPlotter, generate_synthetic_data


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic data and/or plot distributions."
    )
    parser.add_argument(
        "--input", required=True, help="Path to the original dataset (CSV file)."
    )
    parser.add_argument(
        "--output", help="Path to store the generated synthetic data (CSV file)."
    )
    parser.add_argument(
        "--rows", type=int, help="Number of rows for the synthetic data."
    )
    parser.add_argument("--output-dir", help="Directory to save comparison plots.")
    args = parser.parse_args()

    if args.rows and not args.output:
        parser.error("--output must be specified if --rows is specified.")

    # Read the input file
    df = pd.read_csv(args.input, sep="\t")

    synthetic_df = None
    if args.rows:
        # Generate synthetic data
        synthetic_df = generate_synthetic_data(df, args.rows)

        # Save the synthetic data to output file
        if args.output:
            # Ensure output directory exists
            output_dir = os.path.dirname(args.output)
            if not os.path.exists(output_dir) and output_dir:
                os.makedirs(output_dir)
            synthetic_df.to_csv(args.output, index=False)

    # Create an instance of DataPlotter
    plotter = DataPlotter(metadata=df, synthetic_metadata=synthetic_df)

    # Plot numerical distributions
    plotter.plot_numerical_distributors()

    # get a random row from the dataframe
    random_row = df.iloc[np.random.randint(0, df.shape[0])]

    plotter_with_rows = DataPlotter(metadata=df, target_row=random_row)

    # Plot numerical distributions with a target row
    plotter_with_rows.plot_numerical_distributors()

    # Plot comparisons if synthetic data was generated
    if synthetic_df is not None and args.output_dir:
        plotter.plot_comparisons(output_dir=args.output_dir)


if __name__ == "__main__":
    main()
