
# Standard imports
import os
import sys

# Third party imports
import pandas as pd

# Local imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from SOM.utils.SOM import SOM

train_dat = pd.read_csv("../../data/sim_data_pseudo_feature.csv")

other_dat = pd.read_csv("../../data/sim_data_umap_data.csv").iloc[:, [1, 2, 6]]

som = SOM(
    train_dat=train_dat,
    other_dat=other_dat,
    scale_method="zscore",
    x_dim=4,
    y_dim=3,
    topology="hexagonal",
    neighborhood_fnc="gaussian",
    epochs=20
)

som.train_map()

# som.plot_component_planes(
#     output_dir="output_figs"
# )

som.plot_categorical_data(
    output_dir="output_figs"
)
