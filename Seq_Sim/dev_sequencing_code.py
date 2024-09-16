import pandas as pd
import numpy as np
from utils import DataPlotter
import argparse


def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(
        description="Run the genetic sequencing simulation and plot distributions."
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        help="Path to the original genetic sequencing dataset (CSV file).",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Path to store the simulated data (default: None).",
    )
    args = parser.parse_args()

    # Load the data
    df = pd.read_csv(str(args.input), sep="\t")

    # Get a random row from the dataframe
    random_row = df.iloc[np.random.randint(0, df.shape[0])]

    # Create an instance of DataPlotter
    plotter_without_row = DataPlotter(metadata=df)
    plotter_with_row = DataPlotter(metadata=df, target_row=random_row)

    # Plot numerical distributions
    plotter_with_row.plot_numerical_distributors()  # No need to pass target_row here
    plotter_without_row.plot_numerical_distributors()  # This already uses the target_row from initialization

    # Uncomment the following lines to plot categorical distributions
    # plotter.plot_categorical_distributions()
    # plotter.plot_categorical_distributions()  # This already uses the target_row from initialization


if __name__ == "__main__":
    main()
