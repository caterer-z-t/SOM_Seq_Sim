import pandas as pd
import numpy as np
from utils import DataPlotter
import sys


def main():
    input_file = sys.argv[1]

    df = pd.read_csv(str(input_file), sep="\t")

    # Get a random row from the dataframe
    random_row = df.iloc[np.random.randint(0, df.shape[0])]

    # Create an instance of DataPlotter
    plotter = DataPlotter(metadata=df, target_row=random_row)

    # Plot numerical distributions
    plotter.plot_numerical_distributors()  # No need to pass target_row here
    plotter.plot_numerical_distributors()  # This already uses the target_row from initialization

    # Uncomment the following lines to plot categorical distributions
    # plotter.plot_categorical_distributions()
    # plotter.plot_categorical_distributions()  # This already uses the target_row from initialization


if __name__ == "__main__":
    main()
