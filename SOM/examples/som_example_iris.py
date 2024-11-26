"""
Basic Example of Using the Self-Organizing Map (SOM) Class
==========================================================

This script provides a straightforward example of how to use the SOM class to train and evaluate 
a Self-Organizing Map (SOM) with predetermined hyperparameters. The process includes calculating 
fit metrics and visualizing the results.

Key Features:
-------------
1. Predetermined Hyperparameters:
    - The SOM is configured with a fixed set of hyperparameters, including scaling method, grid 
      dimensions, topology, neighborhood function, and number of epochs.

2. Fit Metric Calculation:
    - The Percent Variance Explained (PVE) and Topographic Error are calculated to evaluate the 
      quality of the SOM fit.

3. Visualization:
    - Component planes for individual features are generated to show their distribution across the 
      SOM grid.
    - Categorical data distributions are visualized on the SOM grid to examine relationships 
      between numerical and categorical variables.

Considerations:
---------------
This script serves as a basic example and is not intended for rigorous analysis or hyperparameter 
tuning. For more complex workflows, such as automated hyperparameter optimization, refer to scripts 
specifically designed for grid search or other tuning techniques.
"""

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
