import pytest
import pandas as pd
from SOM_Seq_Sim.Seq_Sim.utils import SimSeq, load_config
import numpy as np
from unittest import mock
import os


def test_optimal_layout():
    # Create an instance of DataPlotter with dummy data
    dummy_metadata = pd.DataFrame()
    plotter = SimSeq(metadata=dummy_metadata)

    assert plotter.optimal_subplot_layout(9) == (3, 3)
    assert plotter.optimal_subplot_layout(10) == (2, 5)
    assert plotter.optimal_subplot_layout(11) == (6, 2)
    assert plotter.optimal_subplot_layout(12) == (2, 6)


def test_remainder():
    # Create an instance of DataPlotter with dummy data
    dummy_metadata = pd.DataFrame()
    plotter = SimSeq(metadata=dummy_metadata)

    assert plotter.remainder(10, 3) == 1
    assert plotter.remainder(12, 5) == 2
    assert plotter.remainder(7, 2) == 1
    assert plotter.remainder(6, 3) == 0


# Mock the plt.savefig to prevent actual file creation during tests
@mock.patch("matplotlib.pyplot.savefig")
def test_plot_numerical_distributors(mock_savefig):
    # Create a DataFrame with some numerical data
    data = pd.DataFrame({"col1": np.random.rand(100), "col2": np.random.rand(100)})

    # Instantiate the DataPlotter
    plotter = SimSeq(metadata=data)

    # Test plotting numerical distributions
    plotter.plot_numerical_distributors(output_dir="./")

    # Ensure the savefig was called to generate the plot
    assert mock_savefig.called


@mock.patch("matplotlib.pyplot.savefig")
def test_plot_categorical_distributors(mock_savefig):
    # Create a DataFrame with some categorical data
    data = pd.DataFrame(
        {
            "cat1": np.random.choice(["A", "B", "C"], size=100),
            "cat2": np.random.choice(["X", "Y", "Z"], size=100),
        }
    )

    # Instantiate the DataPlotter
    plotter = SimSeq(metadata=data)

    # Test plotting categorical distributions
    plotter.plot_categorical_distributions(output_dir="./")

    # Ensure the savefig was called to generate the plot
    assert mock_savefig.called


def test_generate_synthetic_data():
    # Create a DataFrame with both categorical and numerical data
    original_df = pd.DataFrame(
        {
            "numerical": np.random.randn(100),
            "categorical": np.random.choice(["A", "B", "C"], size=100),
        }
    )

    # Instantiate the DataPlotter
    plotter = SimSeq(metadata=original_df)

    # Generate synthetic data
    synthetic_df = plotter.generate_synthetic_data(num_rows=100)

    # Check that the synthetic data has the same columns
    assert list(synthetic_df.columns) == list(original_df.columns)

    # Check that the synthetic data has the correct number of rows
    assert len(synthetic_df) == 100

    # Check that the data types match
    assert synthetic_df["categorical"].dtype == original_df["categorical"].dtype


def test_load_config_from_args():
    # Mock command-line arguments
    args = mock.Mock()
    args.input = "/path/to/input.csv"
    args.output = "/path/to/output.csv"
    args.rows = None
    args.config = None
    args.output_dir = "/path/to/output_dir"

    input_path, output_path, rows, output_dir = load_config(args)

    assert input_path == "/path/to/input.csv"
    assert output_path == "/path/to/output.csv"
    assert rows is None
    assert output_dir == "/path/to/output_dir"


def test_load_config_from_yaml(tmpdir):
    # Create a temporary YAML file for the test
    config_path = tmpdir.join("config.yml")
    config_data = """
    input: /path/to/input.csv
    output: /path/to/output.csv
    rows: 100
    output_dir: /path/to/output_dir
    """
    config_path.write(config_data)

    args = mock.Mock()
    args.config = str(config_path)
    args.input = None
    args.output = None
    args.rows = None
    args.output_dir = None

    input_path, output_path, rows, output_dir = load_config(args)

    assert input_path == "/path/to/input.csv"
    assert output_path == "/path/to/output.csv"
    assert rows == 100
    assert output_dir == "/path/to/output_dir"
