"""
Hyperparameter Tuning Example for Self-Organizing Maps (SOM)
===========================================================

This script demonstrates how to use the SOM class to perform hyperparameter tuning. By defining a
grid of hyperparameters and systematically testing all combinations, the script identifies the
optimal configuration of hyperparameters for a given dataset.

Features of This Script:
-------------------------
1. Grid Search Implementation:
    - Hyperparameters such as `scale_method`, `x_dim`, `y_dim`, `topology`, `neighborhood_fnc`,
      and `epochs` are systematically tested over a predefined range of values.
    - The metrics Percent Variance Explained (PVE) and Topographic Error are combined into
      a scoring function to evaluate the SOM's performance for each combination.

2. Ease of Use:
    - The script leverages Python's `itertools.product` for a clean and systematic exploration of
      hyperparameter combinations.
    - Metrics are calculated using the SOM class' built-in methods, making the evaluation process
      seamless.

3. Visualization of Results:
    - Once the best hyperparameters are identified, the SOM is retrained, and component planes
      and categorical data distributions are visualized and saved.

Considerations:
---------------
- Time Complexity:
    - Depending on the size of the dataset and the number of hyperparameter combinations, this 
      process may take significant time. Adjust the ranges of the hyperparameters to balance 
      between thoroughness and computational efficiency.
    
- Extensibility:
    - The scoring function can be adjusted based on specific requirements. In this example, the 
      score is computed as PVE minus a scaled Topographic Error.

Output:
-------
- Best hyperparameters and their resulting score
- Final SOM trained with the best parameters.
- Saved visualizations of component planes and categorical data distributions.

Usage:
------
Modify the dataset path and hyperparameter grid as needed.
"""

# Standard imports
import itertools
import os
import sys

# Third party imports
import pandas as pd

# Local imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from SOM.utils.som_utils import SOM

# Load data
train_dat = pd.read_csv("../../data/titanic_training_data.csv")
other_dat = pd.read_csv("../../data/titanic_categorical_data.csv")

# Define hyperparameter grid
hyperparameter_grid = {
    "scale_method": ["zscore", "minmax"],
    "x_dim": [3, 5, 7],
    "y_dim": [2, 4, 6],
    "topology": ["rectangular", "hexagonal"],
    "neighborhood_fnc": ["gaussian", "bubble"],
    "epochs": [50, 100, 200],
}

# Initialize variables to store the best parameters and score
best_params = None
best_score = -float("inf")

# Perform grid search
for params in itertools.product(
    hyperparameter_grid["scale_method"],
    hyperparameter_grid["x_dim"],
    hyperparameter_grid["y_dim"],
    hyperparameter_grid["topology"],
    hyperparameter_grid["neighborhood_fnc"],
    hyperparameter_grid["epochs"],
):
    scale_method, x_dim, y_dim, topology, neighborhood_fnc, epochs = params

    # Train SOM with the current hyperparameter combination
    som = SOM(
        train_dat=train_dat,
        other_dat=other_dat,
        scale_method=scale_method,
        x_dim=x_dim,
        y_dim=y_dim,
        topology=topology,
        neighborhood_fnc=neighborhood_fnc,
        epochs=epochs
    )
    som.train_map()

    # Calculate evaluation metrics
    pve = som.calculate_percent_variance_explained()
    topographic_error = som.calculate_topographic_error()

    # Combine metrics into a single score (higher PVE and lower error are better)
    score = pve - topographic_error * 100

    print(f"Tested params: {params} | Score: {score:.2f}")

    # Update the best parameters if the current score is better
    if score > best_score:
        best_score = score
        best_params = params


# Output the best hyperparameters
print("\nBest Hyperparameters:")
print(f"Scale Method: {best_params[0]}")
print(f"x_dim: {best_params[1]}")
print(f"y_dim: {best_params[2]}")
print(f"Topology: {best_params[3]}")
print(f"Neighborhood Function: {best_params[4]}")
print(f"Epochs: {best_params[5]}")
print(f"Best Score: {best_score:.2f}")

# Train and visualize the best SOM
best_som = SOM(
    train_dat=train_dat,
    other_dat=other_dat,
    scale_method=best_params[0],
    x_dim=best_params[1],
    y_dim=best_params[2],
    topology=best_params[3],
    neighborhood_fnc=best_params[4],
    epochs=best_params[5]
)
best_som.train_map()

# Get fit metrics
pve = best_som.calculate_percent_variance_explained()
topographic_error = best_som.calculate_topographic_error()

print("\nFinal SOM Performance:")
print(f"Percent variance explained = {pve}%")
print(f"Topographic error = {topographic_error}")

# Plot component planes
best_som.plot_component_planes(
    output_dir="output/titanic"
)

# Plot SOM Map Using Categorical Data
best_som.plot_categorical_data(
    output_dir="output/titanic"
)
