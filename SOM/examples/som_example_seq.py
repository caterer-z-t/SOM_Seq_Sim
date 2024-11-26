"""
Application of Self-Organizing Maps (SOM) for Cell Type Identification
======================================================================

This script demonstrates the original application for which the SOM class and associated project 
were developed: identifying different cell types based on simulated single-cell data. By training 
a Self-Organizing Map (SOM) on high-dimensional data, this script allows for efficient clustering 
and visualization of cellular features.

Key Features:
-------------
1. Application-Specific Data:
    - The training data represents pseudo-features derived from simulated single-cell experiments.
    - Additional data (`other_dat`) includes the cell type (and other categorical data) labels for
      visualizing categorical clustering results.

2. Visualization of Results:
    - The SOM grid is used to map and visualize cell-type distributions.
    - By examining the output figures, it becomes evident that the SOM successfully identifies 
      distinct cell types based on the input data.

3. Efficient Handling of High-Dimensional Data:
    - The SOM can handle 1000 features in the training data, making it suitable for large-scale 
      single-cell datasets.

Considerations:
---------------
1. Output Figures:
    - Generating component planes for all 1000 features produces a large number of figures. 
      Uncomment the relevant section of the script to enable this visualization if desired.

2. Cell Type Identification:
    - The SOM's ability to separate and group different cell types can be verified by analyzing the 
      categorical data visualizations. Each cell type forms distinct clusters, demonstrating the
      SOM's effectiveness in preserving topology and capturing meaningful patterns in high-
      dimensional data.
"""

# Standard imports
import os
import sys

# Third party imports
import pandas as pd

# Local imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from SOM.utils.som_utils import SOM


# Load data
train_dat = pd.read_csv("../../data/sim_data_pseudo_feature.csv")
other_dat = pd.read_csv("../../data/sim_data_umap_data.csv")

# Train SOM
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

# Train SOM Map
som.train_map()

# Since there are 1000  features, this creates 1000 figures. Uncomment if you want to create
# these figures
# som.plot_component_planes(
#     output_dir="output/seq"
# )

# Plot SOM Map Using Categorical Data
som.plot_categorical_data(
    output_dir="output/seq"
)
