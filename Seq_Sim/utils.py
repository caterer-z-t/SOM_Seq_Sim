import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
from matplotlib.backends.backend_pdf import PdfPages
import math


def generate_synthetic_data(original_df, num_rows):
    """
    Generate a synthetic DataFrame with the same column distributions as the original DataFrame.

    Args:
        original_df (pd.DataFrame): The original DataFrame to base the synthetic data on.
        num_rows (int): The number of rows for the synthetic DataFrame.

    Returns:
        pd.DataFrame: A synthetic DataFrame with the specified number of rows.
    """
    synthetic_df = pd.DataFrame()

    for column in original_df.columns:
        if (
            original_df[column].dtype == "object"
            or original_df[column].dtype.name == "category"
        ):
            # Categorical data
            proportions = original_df[column].value_counts(normalize=True)
            synthetic_df[column] = np.random.choice(
                proportions.index, size=num_rows, p=proportions.values
            )
        else:
            # Numerical data
            kde = gaussian_kde(original_df[column].dropna())
            synthetic_df[column] = kde.resample(num_rows).T.flatten()

    return synthetic_df


class DataPlotter:
    def __init__(self, metadata, target_row=None, synthetic_metadata=None):
        """
        Initializes the DataPlotter with metadata and an optional target row.

        Args:
            metadata (pd.DataFrame): The dataframe containing the metadata.
            target_row (pd.Series, optional): The row containing the values to highlight. Defaults to None.
        """
        self.metadata = metadata
        self.target_row = target_row
        self.synthetic_metadata = synthetic_metadata

    @staticmethod
    def remainder(dividend, divisor):
        """
        Returns the remainder after dividing the dividend by the divisor.

        Args:
            dividend (int or float): The number to be divided.
            divisor (int or float): The number by which to divide.

        Returns:
            int or float: The remainder after division.
        """
        return dividend % divisor

    @staticmethod
    def optimal_subplot_layout(num_plots, max_cols=6):
        """
        Calculate optimal number of rows and columns for subplots to minimize white space.

        Args:
            num_plots (int): Number of plots to display.
            max_cols (int): Maximum number of columns per row.

        Returns:
            tuple: (num_rows, num_cols)
        """
        best_rows = None
        best_cols = None
        min_remainder = float("inf")

        for cols in range(1, max_cols + 1):
            rows = math.ceil(num_plots / cols)
            left_over = DataPlotter.remainder(rows * cols, num_plots)

            if left_over < min_remainder or (
                left_over == min_remainder and cols > best_cols
            ):
                min_remainder = left_over
                best_rows, best_cols = rows, cols

        if best_cols == 1:
            for cols in range(2, max_cols + 1):
                rows = math.ceil(num_plots / cols)
                left_over = DataPlotter.remainder(rows * cols, num_plots)

                if left_over == 1:
                    best_rows, best_cols = rows, cols
                    break

        return best_rows, best_cols

    def plot_numerical_distributors(
        self, output_dir=None, pdf=None, max_plots_per_page=12
    ):
        """
        Plots histograms and KDE with peaks for numerical columns in the metadata.
        Highlights the target row values if provided.

        Args:
            output_dir (str, optional): The directory to save the plot. Defaults to None.
            pdf (PdfPages, optional): PdfPages object to save the plot in a PDF. Defaults to None.
            max_plots_per_page (int, optional): Maximum number of plots to display on each page. Defaults to 12.
        """
        numerical_metadata = self.metadata.select_dtypes(include=["int64", "float64"])

        if numerical_metadata.empty:
            raise ValueError("No numerical columns found in the DataFrame.")

        num_numerical = len(numerical_metadata.columns)
        num_plots = num_numerical * 2

        page_count = (num_plots + max_plots_per_page - 1) // max_plots_per_page

        for page in range(page_count):
            start_idx = page * max_plots_per_page // 2
            end_idx = min(start_idx + max_plots_per_page // 2, num_numerical)

            num_rows, num_cols = self.optimal_subplot_layout((end_idx - start_idx) * 2)

            fig, axs = plt.subplots(num_rows, num_cols, figsize=(18, num_rows * 5))
            axs = axs.flatten()

            for i, column in enumerate(numerical_metadata.columns[start_idx:end_idx]):
                sns.histplot(
                    numerical_metadata[column].dropna(), kde=True, ax=axs[2 * i]
                )

                if self.target_row is not None:
                    target_value = self.target_row[column]
                    axs[2 * i].axvline(
                        target_value,
                        color="red",
                        linestyle="--",
                        label=f"Target: {target_value}",
                    )
                    axs[2 * i].legend()

                axs[2 * i].set_title(f"Histogram & KDE: {column}")

                kde = gaussian_kde(numerical_metadata[column].dropna())
                x = np.linspace(
                    numerical_metadata[column].min(),
                    numerical_metadata[column].max(),
                    1000,
                )
                y = kde(x)
                peaks_kde, _ = find_peaks(y)

                axs[2 * i + 1].plot(x, y, label="KDE", color="blue")
                axs[2 * i + 1].hist(
                    numerical_metadata[column].dropna(),
                    bins=20,
                    density=True,
                    alpha=0.5,
                    color="gray",
                )
                axs[2 * i + 1].plot(x[peaks_kde], y[peaks_kde], "ro", label="Peaks")

                if self.target_row is not None:
                    axs[2 * i + 1].axvline(
                        target_value,
                        color="red",
                        linestyle="--",
                        label=f"Target: {target_value}",
                    )
                    axs[2 * i + 1].legend()

                axs[2 * i + 1].set_title(f"KDE with Peaks: {column}")

            for j in range((end_idx - start_idx) * 2, len(axs)):
                axs[j].axis("off")

            plt.tight_layout(pad=3.0)

            if output_dir:
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                plt.savefig(
                    os.path.join(
                        output_dir, f"numerical_distributions_page_{page + 1}.png"
                    )
                )
            elif pdf:
                pdf.savefig()
            else:
                plt.show()

            plt.close()

    def plot_categorical_distributions(
        self, output_dir=None, pdf=None, max_plots_per_page=12
    ):
        """
        Plots bar charts and pie charts for categorical columns in the metadata.
        Highlights the target row values if provided.

        Args:
            output_dir (str, optional): The directory to save the plot. Defaults to None.
            pdf (PdfPages, optional): PdfPages object to save the plot in a PDF. Defaults to None.
            max_plots_per_page (int, optional): Maximum number of plots to display on each page. Defaults to 12.
        """
        categorical_metadata = self.metadata.select_dtypes(
            include=["object", "category"]
        )

        if categorical_metadata.empty:
            raise ValueError("No categorical columns found in the DataFrame.")

        num_categorical = len(categorical_metadata.columns)
        num_plots = num_categorical * 2

        page_count = (num_plots + max_plots_per_page - 1) // max_plots_per_page

        for page in range(page_count):
            start_idx = page * max_plots_per_page // 2
            end_idx = min(start_idx + max_plots_per_page // 2, num_categorical)

            num_rows, num_cols = self.optimal_subplot_layout((end_idx - start_idx) * 2)

            fig, axs = plt.subplots(num_rows, num_cols, figsize=(18, num_rows * 5))
            axs = axs.flatten()

            for i, column in enumerate(categorical_metadata.columns[start_idx:end_idx]):
                proportions = self.metadata[column].value_counts(normalize=True)

                sns.barplot(
                    x=proportions.index,
                    y=proportions.values,
                    ax=axs[2 * i],
                    palette="viridis",
                )
                axs[2 * i].set_ylabel("Proportion")
                axs[2 * i].set_xticklabels(
                    axs[2 * i].get_xticklabels(), rotation=45, ha="right"
                )

                axs[2 * i + 1].pie(
                    proportions,
                    labels=proportions.index,
                    autopct="%1.1f%%",
                    colors=sns.color_palette("viridis", len(proportions)),
                )

                if self.target_row is not None:
                    target_value = self.target_row[column]
                    if target_value in proportions.index:
                        target_idx = list(proportions.index).index(target_value)
                        axs[2 * i].axvline(
                            target_idx,
                            color="red",
                            linestyle="--",
                            label=f"Target: {target_value}",
                        )
                        axs[2 * i].legend()

            for j in range((end_idx - start_idx) * 2, len(axs)):
                axs[j].axis("off")

            plt.tight_layout(pad=3.0)

            if output_dir:
                if not os.path.exists(output_dir):
                    os.makedirs(output_dir)
                plt.savefig(
                    os.path.join(
                        output_dir, f"categorical_distributions_page_{page + 1}.png"
                    )
                )
            elif pdf:
                pdf.savefig()
            else:
                plt.show()

            plt.close()

    def generate_synthetic_data(self, num_rows):
        """
        Generates a new DataFrame with the same distributions or proportions as the original dataset.

        Args:
            num_rows (int): Number of rows for the new DataFrame.

        Returns:
            pd.DataFrame: A new DataFrame with synthetic data.
        """
        synthetic_data = pd.DataFrame()

        # Numerical columns
        numerical_metadata = self.metadata.select_dtypes(include=["int64", "float64"])
        for column in numerical_metadata.columns:
            data = numerical_metadata[column].dropna()
            kde = gaussian_kde(data)
            x = np.linspace(data.min(), data.max(), 1000)
            y = kde(x)
            synthetic_data[column] = np.random.choice(x, size=num_rows, p=y / y.sum())

        # Categorical columns
        categorical_metadata = self.metadata.select_dtypes(
            include=["object", "category"]
        )
        for column in categorical_metadata.columns:
            proportions = self.metadata[column].value_counts(normalize=True)
            synthetic_data[column] = np.random.choice(
                proportions.index, size=num_rows, p=proportions.values
            )

        return synthetic_data

    def plot_comparisons(self, output_dir=None, pdf=None):
        """
        Plots comparisons between the original and synthetic data for numerical columns.

        Args:
            output_dir (str, optional): The directory to save the plot. Defaults to None.
        """
        if self.synthetic_metadata is None:
            raise ValueError("No synthetic data provided for comparison.")

        numerical_metadata = self.metadata.select_dtypes(include=["int64", "float64"])
        numerical_synthetic = self.synthetic_metadata.select_dtypes(
            include=["int64", "float64"]
        )

        if numerical_metadata.empty or numerical_synthetic.empty:
            raise ValueError("No numerical columns found in one or both DataFrames.")

        num_numerical = len(numerical_metadata.columns)
        num_plots = num_numerical * 2

        num_rows, num_cols = self.optimal_subplot_layout(num_plots)

        fig, axs = plt.subplots(num_rows, num_cols, figsize=(18, num_rows * 5))
        axs = axs.flatten()

        for i, column in enumerate(numerical_metadata.columns):
            # Plot original data
            sns.histplot(
                numerical_metadata[column].dropna(),
                kde=True,
                ax=axs[2 * i],
                color="blue",
                label="Original",
            )
            sns.histplot(
                numerical_synthetic[column].dropna(),
                kde=True,
                ax=axs[2 * i],
                color="red",
                label="Synthetic",
                alpha=0.5,
            )

            axs[2 * i].set_title(f"Histogram & KDE: {column}")
            axs[2 * i].legend()

            # Plot KDE with peaks
            kde_original = gaussian_kde(numerical_metadata[column].dropna())
            kde_synthetic = gaussian_kde(numerical_synthetic[column].dropna())
            x = np.linspace(
                min(
                    numerical_metadata[column].min(), numerical_synthetic[column].min()
                ),
                max(
                    numerical_metadata[column].max(), numerical_synthetic[column].max()
                ),
                1000,
            )
            y_original = kde_original(x)
            y_synthetic = kde_synthetic(x)
            peaks_original, _ = find_peaks(y_original)
            peaks_synthetic, _ = find_peaks(y_synthetic)

            axs[2 * i + 1].plot(x, y_original, label="Original KDE", color="blue")
            axs[2 * i + 1].plot(x, y_synthetic, label="Synthetic KDE", color="red")
            axs[2 * i + 1].hist(
                numerical_metadata[column].dropna(),
                bins=20,
                density=True,
                alpha=0.5,
                color="blue",
            )
            axs[2 * i + 1].hist(
                numerical_synthetic[column].dropna(),
                bins=20,
                density=True,
                alpha=0.5,
                color="red",
            )
            axs[2 * i + 1].plot(
                x[peaks_original],
                y_original[peaks_original],
                "bo",
                label="Original Peaks",
            )
            axs[2 * i + 1].plot(
                x[peaks_synthetic],
                y_synthetic[peaks_synthetic],
                "ro",
                label="Synthetic Peaks",
            )

            axs[2 * i + 1].set_title(f"KDE with Peaks: {column}")
            axs[2 * i + 1].legend()

        for j in range(num_plots, len(axs)):
            axs[j].axis("off")

        plt.tight_layout(pad=3.0)

        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            plt.savefig(os.path.join(output_dir, "comparison_plots.png"))
        elif pdf:
            pdf.savefig()
        else:
            plt.show()

        plt.close()
