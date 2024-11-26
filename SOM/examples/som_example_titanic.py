# Standard imports
import os
import sys

# Third party imports
import pandas as pd

# Local imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from SOM.utils.som import SOM

# Load data
data = pd.read_csv("../../data/titanic.csv")
#data=pd.read_csv("data/titanic.csv")
data.dropna(
    axis=0,
    how='any',
    subset=None,
    inplace=True
)
train_dat = data.iloc[:, 0:5]
other_dat = data.iloc[:, 5:]

# Train SOM
som = SOM(
    train_dat=train_dat,
    other_dat=other_dat,
    scale_method="zscore",
    x_dim=5,
    y_dim=4,
    topology="hexagonal",
    neighborhood_fnc="gaussian",
    epochs=100
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
    output_dir="output_figs/titanic"
)

# Plot SOM Map Using Categorical Data
som.plot_categorical_data(
    output_dir="output_figs/titanic"
)
