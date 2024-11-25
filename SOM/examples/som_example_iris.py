# Standard imports
import os
import sys

# Third party imports
import pandas as pd

# Local imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from SOM.utils.som import SOM

# Load data
data = pd.read_csv("../../data/iris.csv")
#data = pd.read_csv("data/iris.csv")
train_dat = data.iloc[:, 0:4]
other_dat = data.iloc[:, 4].to_frame()

# Train SOM
som = SOM(
    train_dat=train_dat,
    other_dat=other_dat,
    scale_method="minmax",
    x_dim=5,
    y_dim=2,
    topology="hexagonal",
    neighborhood_fnc="gaussian",
    epochs=17
)

# Train SOM Map
som.train_map()

# Get fit metrics
pve = som.calculate_percent_variance_explained()
topographic_error = som.calculate_topographic_error()

print(f"Percent variance explained = {pve}%")
print(f"Topographic error = {topographic_error}")

# Plot component planes
som.plot_component_planes(
    output_dir="output_figs/iris"
)

# Plot SOM Map Using Categorical Data
som.plot_categorical_data(
    output_dir="output_figs/iris"
)
