import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from scipy.signal import find_peaks
from matplotlib.backends.backend_pdf import PdfPages
import math


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


def optimal_subplot_layout(num_plots, max_cols=6):
    """
    Calculate optimal number of rows and columns for subplots to minimize white space.
    Prioritize the layout with the highest number of columns in case of ties in leftover space.
    If the optimal layout has only 1 column, pick the best layout that has 1 remainder with maximum columns.

    Args:
        num_plots (int): Number of plots to display.
        max_cols (int): Maximum number of columns per row.

    Returns:
        tuple: (num_rows, num_cols)
    """
    best_rows = None
    best_cols = None
    min_remainder = float("inf")  # Start with infinity as the minimum remainder

    for cols in range(1, max_cols + 1):
        rows = math.ceil(num_plots / cols)
        left_over = remainder(
            rows * cols, num_plots
        )  # Calculate how much space is left unused

        # Update if this combination produces fewer leftover spaces,
        # or if the remainder is the same but the number of columns is higher.
        if left_over < min_remainder or (
            left_over == min_remainder and cols > best_cols
        ):
            min_remainder = left_over
            best_rows, best_cols = rows, cols

    # Special case: If the best option results in only 1 column, find a better option
    if best_cols == 1:
        for cols in range(2, max_cols + 1):
            rows = math.ceil(num_plots / cols)
            left_over = remainder(rows * cols, num_plots)

            # Prefer layouts that leave exactly 1 remainder, with the most columns
            if left_over == 1:
                best_rows, best_cols = rows, cols
                break  # Stop as soon as we find a layout with 1 remainder

    return best_rows, best_cols


def plot_numerical_distributions(
    metadata, output_dir=None, pdf=None, max_plots_per_page=12
):
    """
    Plots histograms with KDE and KDE with peaks for numerical columns in the provided metadata,
    splitting across multiple pages in the PDF if necessary.

    Args:
        metadata (pd.DataFrame): The dataframe containing the metadata with numerical columns.
        output_dir (str, optional): The directory to save the plot. Defaults to None.
        pdf (PdfPages, optional): PdfPages object to save the plot in a PDF. Defaults to None.
        max_plots_per_page (int, optional): Maximum number of plots to display on each page. Defaults to 12.
    """
    # Ensure the dataframe only contains numerical columns
    numerical_metadata = metadata.select_dtypes(include=["int64", "float64"])

    if numerical_metadata.empty:
        raise ValueError("No numerical columns found in the DataFrame.")

    num_numerical = len(numerical_metadata.columns)
    num_plots = num_numerical * 2  # Two plots per column: histogram + KDE with peaks

    page_count = (num_plots + max_plots_per_page - 1) // max_plots_per_page

    for page in range(page_count):
        # Determine the range of plots to display on the current page
        start_idx = page * max_plots_per_page // 2
        end_idx = min(start_idx + max_plots_per_page // 2, num_numerical)

        num_rows, num_cols = optimal_subplot_layout((end_idx - start_idx) * 2)

        fig, axs = plt.subplots(num_rows, num_cols, figsize=(18, num_rows * 5))
        axs = axs.flatten()

        for i, column in enumerate(numerical_metadata.columns[start_idx:end_idx]):
            # Plot histogram with KDE
            sns.histplot(numerical_metadata[column].dropna(), kde=True, ax=axs[i])
            axs[i].set_title(f"Histogram & KDE: {column}")

            # Calculate KDE and peaks for the second row
            kde = gaussian_kde(numerical_metadata[column].dropna())
            x = np.linspace(
                numerical_metadata[column].min(), numerical_metadata[column].max(), 1000
            )
            y = kde(x)
            peaks_kde, _ = find_peaks(y)

            # Plot KDE with peaks
            axs[i + (end_idx - start_idx)].plot(x, y, label="KDE", color="blue")
            axs[i + (end_idx - start_idx)].hist(
                numerical_metadata[column].dropna(),
                bins=20,
                density=True,
                alpha=0.5,
                color="gray",
                label="Histogram",
            )
            axs[i + (end_idx - start_idx)].plot(
                x[peaks_kde], y[peaks_kde], "ro", label="Peaks"
            )
            axs[i + (end_idx - start_idx)].set_title(f"KDE with Peaks: {column}")
            axs[i + (end_idx - start_idx)].set_xlabel(column)
            axs[i + (end_idx - start_idx)].set_ylabel("Density")
            axs[i + (end_idx - start_idx)].legend()

        # Hide any unused subplots
        for j in range((end_idx - start_idx) * 2, len(axs)):
            axs[j].axis("off")

        plt.tight_layout(pad=3.0)  # Adjust padding between subplots

        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            plt.savefig(
                os.path.join(output_dir, f"numerical_distributions_page_{page + 1}.png")
            )
        elif pdf:
            pdf.savefig()  # Add the figure as a new page to the PDF
        else:
            plt.show()

        plt.close()


def plot_categorical_distributions(
    metadata, output_dir=None, pdf=None, max_plots_per_page=12
):
    """
    Plots bar charts and pie charts (as proportions) for categorical columns in the provided metadata,
    splitting across multiple pages in the PDF if necessary.

    Args:
        metadata (pd.DataFrame): The dataframe containing the metadata with categorical columns.
        output_dir (str, optional): The directory to save the plot. Defaults to None.
        pdf (PdfPages, optional): PdfPages object to save the plot in a PDF. Defaults to None.
        max_plots_per_page (int, optional): Maximum number of plots to display on each page. Defaults to 12.
    """
    # Ensure the dataframe only contains categorical columns
    categorical_metadata = metadata.select_dtypes(include=["object", "category"])

    if categorical_metadata.empty:
        raise ValueError("No categorical columns found in the DataFrame.")

    num_categorical = len(categorical_metadata.columns)
    num_plots = num_categorical * 2  # Two plots per column: bar chart + pie chart

    # Calculate number of pages
    page_count = (num_plots + max_plots_per_page - 1) // max_plots_per_page

    for page in range(page_count):
        # Determine the range of plots to display on the current page
        start_idx = page * max_plots_per_page // 2
        end_idx = min(start_idx + max_plots_per_page // 2, num_categorical)

        num_rows, num_cols = optimal_subplot_layout((end_idx - start_idx) * 2)

        fig, axs = plt.subplots(num_rows, num_cols, figsize=(18, num_rows * 5))
        axs = axs.flatten()

        for i, column in enumerate(categorical_metadata.columns[start_idx:end_idx]):
            # Calculate proportions
            proportions = metadata[column].value_counts(normalize=True)

            # Plot bar chart (proportions)
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

            # Plot pie chart (proportions)
            axs[2 * i + 1].pie(
                proportions,
                labels=proportions.index,
                autopct="%1.1f%%",
                colors=sns.color_palette("viridis", len(proportions)),
            )

        # Hide any unused subplots
        for j in range((end_idx - start_idx) * 2, len(axs)):
            axs[j].axis("off")

        plt.tight_layout(pad=3.0)  # Adjust padding between subplots

        if output_dir:
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)
            plt.savefig(
                os.path.join(
                    output_dir, f"categorical_distributions_page_{page + 1}.png"
                )
            )
        elif pdf:
            pdf.savefig()  # Add the figure as a new page to the PDF
        else:
            plt.show()

        plt.close()
