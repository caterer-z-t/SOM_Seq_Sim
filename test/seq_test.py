import pytest
import pandas as pd
from SOM_Seq_Sim.Seq_Sim.utils import DataPlotter


def test_optimal_layout():
    # Create an instance of DataPlotter with dummy data
    dummy_metadata = pd.DataFrame()
    plotter = DataPlotter(metadata=dummy_metadata)

    assert plotter.optimal_subplot_layout(9) == (3, 3)
    assert plotter.optimal_subplot_layout(10) == (2, 5)
    assert plotter.optimal_subplot_layout(11) == (6, 2)
    assert plotter.optimal_subplot_layout(12) == (2, 6)


def test_remainder():
    # Create an instance of DataPlotter with dummy data
    dummy_metadata = pd.DataFrame()
    plotter = DataPlotter(metadata=dummy_metadata)

    assert plotter.remainder(10, 3) == 1
    assert plotter.remainder(12, 5) == 2
    assert plotter.remainder(7, 2) == 1
    assert plotter.remainder(6, 3) == 0
