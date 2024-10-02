import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde


class TimeSeriesSimSeq:
    def __init__(self, df):
        """
        Initialize the TimeSeriesSimSeq object with a DataFrame.

        Args:
            df (pd.DataFrame): The original dataset.
        """
        self.df = df

    def generate_synthetic_time_series(
        self, num_rows, num_time_points=1, group_col=None
    ):
        """
        Generates synthetic time-series data.

        Args:
            num_rows (int): Number of rows (instances) for the synthetic dataset.
            num_time_points (int): Number of time points for each instance in the time series. Defaults to 1.
            group_col (str): The column name that contains the group labels (e.g., 'group' with 'disease' and 'healthy').

        Returns:
            list of pd.DataFrame: A list of DataFrames where each DataFrame corresponds to a time point.
        """
        # List to store the synthetic DataFrames for each time point
        synthetic_data_list = []

        # Numerical columns
        numerical_df = self.df.select_dtypes(include=["int64", "float64"])
        # Categorical columns
        categorical_df = self.df.select_dtypes(include=["object", "category"])

        # Generate synthetic group labels if group_col is provided
        if group_col and group_col in self.df.columns:
            group_proportions = self.df[group_col].value_counts(normalize=True)
            synthetic_groups = np.random.choice(
                group_proportions.index, size=num_rows, p=group_proportions.values
            )
        else:
            synthetic_groups = np.random.choice(range(num_rows), size=num_rows)

        for t in range(num_time_points):
            # Create a new DataFrame for the current time point
            synthetic_data = pd.DataFrame()

            # Generate synthetic data for numerical columns
            for column in numerical_df.columns:
                data = numerical_df[column].dropna()
                kde = gaussian_kde(data)
                x = np.linspace(data.min(), data.max(), 1000)
                y = kde(x)
                synthetic_values = np.random.choice(x, size=num_rows, p=y / y.sum())
                synthetic_data[column] = synthetic_values

            # Generate synthetic data for categorical columns
            for column in categorical_df.columns:
                proportions = categorical_df[column].value_counts(normalize=True)
                synthetic_data[column] = np.random.choice(
                    proportions.index, size=num_rows, p=proportions.values
                )

            # Add group labels to the synthetic data for this time point
            if group_col:
                synthetic_data[group_col] = synthetic_groups

            # Append the DataFrame for this time point to the list
            synthetic_data_list.append(synthetic_data)

        return synthetic_data_list
