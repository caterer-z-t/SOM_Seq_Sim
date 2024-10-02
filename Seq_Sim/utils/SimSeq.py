###############################################################
###
###                     Imports
###
###############################################################


import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
import math


###############################################################
###
###                     Object
###
###############################################################


class SimSeq:
    def __init__(self, df, target_row=None, synthetic_df=None):
        """Initialize the SimSeq object with metadata, target row, and synthetic metadata."""
        try:
            if not isinstance(df, pd.DataFrame):
                raise TypeError("The input df should be a pandas DataFrame.")
            if synthetic_df is not None and not isinstance(synthetic_df, pd.DataFrame):
                raise TypeError("The synthetic_df should be a pandas DataFrame.")

            self.df = df
            self.target_row = target_row
            self.synthetic_df = synthetic_df
        except Exception as e:
            print(f"Error initializing SimSeq: {e}")

    @staticmethod
    def remainder(dividend, divisor):
        try:
            return dividend % divisor
        except ZeroDivisionError:
            print("Divisor cannot be zero.")
            return None
        except Exception as e:
            print(f"An error occurred during remainder calculation: {e}")
            return None

    @staticmethod
    def optimal_subplot_layout(num_plots, max_cols=6):
        try:
            best_rows, best_cols, min_remainder = None, None, float("inf")
            for cols in range(1, max_cols + 1):
                rows = math.ceil(num_plots / cols)
                left_over = SimSeq.remainder(rows * cols, num_plots)
                if left_over is not None and (
                    left_over < min_remainder
                    or (left_over == min_remainder and cols > best_cols)
                ):
                    min_remainder, best_rows, best_cols = left_over, rows, cols

            return best_rows, best_cols
        except Exception as e:
            print(f"Error in subplot layout calculation: {e}")
            return 1, 1

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
        try:
            numerical_df = self.df.select_dtypes(include=["int64", "float64"])

            if numerical_df.empty:
                raise ValueError("No numerical columns found in the DataFrame.")

            num_numerical = len(numerical_df.columns)
            num_plots = num_numerical * 2

            page_count = (num_plots + max_plots_per_page - 1) // max_plots_per_page

            for page in range(page_count):
                start_idx = page * max_plots_per_page // 2
                end_idx = min(start_idx + max_plots_per_page // 2, num_numerical)

                num_rows, num_cols = self.optimal_subplot_layout(
                    (end_idx - start_idx) * 2
                )

                fig, axs = plt.subplots(num_rows, num_cols, figsize=(18, num_rows * 5))
                axs = axs.flatten()

                for i, column in enumerate(numerical_df.columns[start_idx:end_idx]):
                    sns.histplot(numerical_df[column].dropna(), kde=True, ax=axs[2 * i])

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

                    kde = gaussian_kde(numerical_df[column].dropna())
                    x = np.linspace(
                        numerical_df[column].min(),
                        numerical_df[column].max(),
                        1000,
                    )
                    y = kde(x)
                    peaks_kde, _ = find_peaks(y)

                    axs[2 * i + 1].plot(x, y, label="KDE", color="blue")
                    axs[2 * i + 1].hist(
                        numerical_df[column].dropna(),
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
        except Exception as e:
            print(f"Error plotting numerical distributions: {e}")

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
        try:
            categorical_df = self.df.select_dtypes(include=["object", "category"])

            if categorical_df.empty:
                raise ValueError("No categorical columns found in the DataFrame.")

            num_categorical = len(categorical_df.columns)
            num_plots = num_categorical * 2

            page_count = (num_plots + max_plots_per_page - 1) // max_plots_per_page

            for page in range(page_count):
                start_idx = page * max_plots_per_page // 2
                end_idx = min(start_idx + max_plots_per_page // 2, num_categorical)

                num_rows, num_cols = self.optimal_subplot_layout(
                    (end_idx - start_idx) * 2
                )

                fig, axs = plt.subplots(num_rows, num_cols, figsize=(18, num_rows * 5))
                axs = axs.flatten()

                for i, column in enumerate(categorical_df.columns[start_idx:end_idx]):
                    proportions = self.df[column].value_counts(normalize=True)

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
        except Exception as e:
            print(f"Error plotting categorical distributions: {e}")

    def plot_comparisons(self, output_dir=None, pdf=None):
        """
        Plots comparisons between the original and synthetic data for numerical columns.

        Args:
            output_dir (str, optional): The directory to save the plot. Defaults to None.
            pdf (PdfPages, optional): PdfPages object to save the plot in a PDF. Defaults to None.
        """
        try:
            if self.synthetic_df is None:
                raise ValueError("No synthetic data provided for comparison.")

            numerical_df = self.df.select_dtypes(include=["int64", "float64"])
            numerical_synthetic = self.synthetic_df.select_dtypes(
                include=["int64", "float64"]
            )

            if numerical_df.empty or numerical_synthetic.empty:
                raise ValueError(
                    "No numerical columns found in one or both DataFrames."
                )

            num_numerical = len(numerical_df.columns)
            num_plots = num_numerical * 2

            num_rows, num_cols = self.optimal_subplot_layout(num_plots)

            fig, axs = plt.subplots(num_rows, num_cols, figsize=(18, num_rows * 5))
            axs = axs.flatten()

            for i, column in enumerate(numerical_df.columns):
                # Plot original data
                sns.histplot(
                    numerical_df[column].dropna(),
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
                kde_original = gaussian_kde(numerical_df[column].dropna())
                kde_synthetic = gaussian_kde(numerical_synthetic[column].dropna())
                x = np.linspace(
                    min(numerical_df[column].min(), numerical_synthetic[column].min()),
                    max(numerical_df[column].max(), numerical_synthetic[column].max()),
                    1000,
                )
                y_original = kde_original(x)
                y_synthetic = kde_synthetic(x)
                peaks_original, _ = find_peaks(y_original)
                peaks_synthetic, _ = find_peaks(y_synthetic)

                axs[2 * i + 1].plot(x, y_original, label="Original KDE", color="blue")
                axs[2 * i + 1].plot(x, y_synthetic, label="Synthetic KDE", color="red")
                axs[2 * i + 1].hist(
                    numerical_df[column].dropna(),
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
        except Exception as e:
            print(f"Error plotting comparisons: {e}")
