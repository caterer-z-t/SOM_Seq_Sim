"""
Self-Organizing Maps (SOMs) Module
===================================

What is a Self-Organizing Map?
------------------------------
A Self-Organizing Map (SOM) is a type of unsupervised artificial neural network used to group and 
visualize data points that are similar to each other. What makes SOMs unique is that they visualize
these groups (or clusters) in simple 2D grid configuration, even though the original data that they
represent is high-dimensional. This makes SOMs very useful for exploring and discovering hidden
patterns with complex datasets.

How Does a SOM Work?
--------------------
- The 2D Grid: A SOM is visualized as a grid of "neurons," arranged in rows and columns. Each
  neuron represents a part of the dataset being modeled.
  
- Weight Vectors: Behind the scenes, each neuron is associated with a "weight vector." This weight
  vector has the same number of dimensions as the original data. Think of the weight vector as the
  neuron's "position" in the high-dimensional space, summarizing the features of the data it
  represents.

- Training the SOM: During training, the SOM adjusts the weight vectors of the neurons to better
  represent the dataset. When a neuron's weight vector is updated, its neighboring neurons on the
  grid are also updated. This ensures that similar data points are mapped to nearby neurons on the
  grid, preserving the topology of the dataset.

- Mapping Data Points to the Grid: For each data point, the SOM calculates how "close" it is to
  each neuron's weight vector using a distance metric (e.g., Euclidean distance). The data point is
  then assigned to the neuron with the closest weight vector. This process transforms the high-
  dimensional data points into positions on the 2D grid.

How Does a SOM Compare to Other Clustering Methods?
---------------------------------------------------
SOMs are similar to k-means clustering in that they group data points by assigning them to a central
representative (a "neuron" in SOMs, or a "centroid" in k-means). However, SOMs are more structured
because neurons are part of a grid, and their updates are influenced by neighboring neurons. This
constraint creates a smooth, organized map where similar clusters are close to each other, making
SOMs ideal for data visualization.
"""

import os
from typing import Dict, List, Tuple, Union

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import matplotlib.figure
from minisom import MiniSom
import numpy as np
import pandas as pd
import seaborn as sns


class SOM():
    """
    A class for building, training, and visualizing Self-Organizing Maps (SOMs). This class
    simplifies the creation and usage of SOMs, providing configurable parameters for training and
    tools for data scaling and visualization of the trained SOM.

    Attributes (set before training):
    ---------------------------------
     - train_dat (np.ndarray): The raw (unscaled) training data used to fit the SOM.
     - other_dat (np.ndarray): Additional, non-numeric data associated with each data point in 
      `train_dat`, but not used to train the SOM.
     - xdim (int): The number of neurons in the x-dimension of the SOM grid.
     - ydim (int): The number of neurons in the y-dimension of the SOM grid.
     - topology (str): The topology of the SOM grid ('rectangular' or 'hexagonal', associated
       with 4 and 6 neighbors, respectively).
     - neighborhood_fnc (str): The neighborhood function used during training.
     - epochs (int): The number of training epochs (i.e., how many times the dataset is presented 
       during the training process).

    Attributes (set after training):
    --------------------------------
     - map (MiniSom): The trained MiniSom instance, created after calling `train_map`.
     - observation_mapping (np.ndarray): The mapping of data points to their closest neurons.
     - neuron_coordinates (pd.DataFrame): The coordinates of neurons in the SOM grid.
     - weights (pd.DataFrame): The unscaled weights of each neuron after training.
     - weights_scaled (pd.DataFrame): The scaled weights of each neuron after training.
"""

    def __init__(
        self,
        train_dat: Union[pd.DataFrame, np.ndarray],
        other_dat: Union[pd.DataFrame, np.ndarray],
        scale_method: str,
        x_dim: int,
        y_dim: int,
        topology: str,
        neighborhood_fnc: str,
        epochs: int
    ):
        # Convert input data to numpy array
        self.train_dat = train_dat.to_numpy() if isinstance(train_dat, pd.DataFrame) else train_dat
        self.other_dat = other_dat.to_numpy() if isinstance(other_dat, pd.DataFrame) else other_dat

        # Set scaling and unscaling methods
        self._scale = getattr(self, f"_{scale_method}_scale", None)
        self._unscale = getattr(self, f"_{scale_method}_unscale", None)

        self.train_dat_scaled, self._scaling_factors = self._scale()
        self.xdim = x_dim
        self.ydim = y_dim
        self.topology = topology
        self.neighborhood_fnc = neighborhood_fnc
        self.epochs = epochs

        # Attributes to be set later
        self.map = None
        self.observation_mapping = None
        self.neuron_coordinates = None
        self.weights = None
        self.weights_scaled = None


    def train_map(self):
        """
        Train the Self-Organizing Map on the scaled training data.

        This method initializes the SOM using PCA-based weight initialization and trains it using
        the specified parameters. It also calculates neuron weights and maps observations to their
        corresponding neurons.

        Sets the following attributes:
            - `map`: The trained MiniSom instance.
            - `observation_mapping`: Mapping of observations to neurons.
            - `neuron_coordinates`: Coordinates of neurons in the SOM grid.
            - `weights_scaled`: Scaled weights of each neuron.
            - `weights`: Unscaled weights of each neuron.
        """

        n_samples = self.train_dat_scaled.shape[0]
        n_features = self.train_dat_scaled.shape[1]

        som = MiniSom(
            x=self.xdim,  # number of neurons in x-dimension
            y=self.ydim,  # number of neurons in y-dimension
            input_len=n_features,  # number of features/variables describing each observation
            sigma=1, #related to the how many neighbors are considered during fitting process
            learning_rate=0.5,  # define how much weight is adjusted
            neighborhood_function=self.neighborhood_fnc,  # form of the neighbordhood function
            topology=self.topology,  # rectangular or hexagonal (4 or 6 neighbors, respectively)
            activation_distance='euclidean',  # method for distance calculation
            random_seed=0
        )

        # Principle_component_analysis for weight initialization
        som.pca_weights_init(self.train_dat_scaled)

        # SOM training
        som.train(
            data=self.train_dat_scaled,
            num_iteration=self.epochs * n_samples,
            random_order=True,
            verbose=False
        )

        # Set attributes yielded from training
        self.map = som
        self.observation_mapping = self._get_observation_neuron_mappings()
        self.neuron_coordinates = self._get_neuron_coordinates()
        self.weights_scaled = self._get_weights()
        self.weights = self._unscale(self.weights_scaled)


    def _zscore_scale(self):
        """
        Scale data using z-score normalization.

        Returns:
            Tuple[np.ndarray, List[np.ndarray]]:
                - (np.ndarray): The scaled data
                - (List[np.ndarray]): A list containing the means and standard deviations used for
                  scaling.
        """

        means = np.mean(self.train_dat, axis=0)
        stds = np.std(self.train_dat, axis=0)
        return (self.train_dat - means) / stds, [means, stds]


    def _minmax_scale(self):
        """
        Scale data using min/max normalization.

        Returns:
            Tuple[np.ndarray, List[np.ndarray]]:
                - np.ndarray: The scaled data
                - List[np.ndarray]: A list containing the minimum and maximum values used for
                  scaling.
        """

        min_vals = np.min(self.train_dat, axis=0)
        max_vals = np.max(self.train_dat, axis=0)
        return (self.train_dat - min_vals) / (max_vals - min_vals), [min_vals, max_vals]


    def _zscore_unscale(self, scaled_data):
        """
        Unscale z-score normalized data back to its original scale.

        Args:
            scaled_data (np.ndarray): Scaled data, often neuron weights returned from the training
            process.

        Returns:
            np.ndarray: The unscaled data.
        """

        means, stds = self._scaling_factors
        return scaled_data * stds + means


    def _minmax_unscale(self, scaled_data):
        """
        Unscale z-score normalized data back to its original scale.

        Args:
            scaled_data (np.ndarray): Scaled data, often neuron weights returned from the training
            process.

        Returns:
            np.ndarray: The unscaled data.
        """

        min_vals, max_vals = self._scaling_factors
        return scaled_data * (max_vals - min_vals) + min_vals


    def _get_observation_neuron_mappings(self):
        """
        Get the neuron id that each observation in `training_dat` is assigned to.

        Returns:
            np.ndarray: An array where each element corresponds to the ID of the winning neuron for
            the respective data point.
        """

        # Get the x-y coordiantes of the winning neurons
        winner_coordinates = np.array([
            self.map.winner(x) for x in self.train_dat_scaled
        ]).T

        # Convert x-y coordinate to "neuron id"
        winner_neuron = np.ravel_multi_index(
            multi_index=winner_coordinates,
            dims=(self.xdim, self.ydim)
        )
        return winner_neuron


    def _get_neuron_coordinates(self):
        """
        Retrieve the coordinates of all neurons in the SOM grid.

        Returns:
            pd.DataFrame: A DataFrame containing the x and y coordinates of each neuron, where the
            first column corresponds to x and the second column corresponds to y.
        """

        HEX_CORRECTION = np.sqrt(3)/2

        xx, yy = self.map.get_euclidean_coordinates()

        x_coords = [xx[(i, j)] for i in range(self.xdim) for j in range(self.ydim)]
        if self.topology == 'hexagonal':
            y_coords = [
                yy[(i, j)] * HEX_CORRECTION for i in range(self.xdim) for j in range(self.ydim)
            ]
        else:
            y_coords = [yy[(i, j)] for i in range(self.xdim) for j in range(self.ydim)]

        return pd.DataFrame({
            'x': x_coords,
            'y': y_coords
        })


    def _get_weights(self):
        """
        Retrieve the weights of each neuron after training.

        Returns:
            pd.DataFrame: A DataFrame where each row corresponds to a neuron's weight vector and 
            the number of columns equal the dimensionality of the original dataset.
        """

        weights_scaled = [
            self.map.get_weights()[i, j, :] for i in range(self.xdim) for j in range(self.ydim)
        ]
        return pd.DataFrame(weights_scaled)


    def plot_component_planes(
        self,
        output_dir: str
    ):
        """
        Generate and save component plane plots for each feature in the training data.

        A component plane is a visualization that shows how the values of a specific feature 
        in `training_dat` are distributed across the SOM grid. Each neuron on the grid is assigned
        a color based on the value of the feature it represents, allowing you to see patterns and
        relationships in the data.

        Component planes help in understanding how individual features contribute to the
        organization of the SOM. By examining these plots, you can:
            - Identify clusters and trends in specific features.
            - Compare how different features vary across the grid, potentially highlighting
              relationships/correlations between features.

        Args:
            - output_dir (str): The directory where the component plane plots will be saved. 

        Saves:
            - One plot per feature as a .png file in the specified directory. Each plot displays
              the SOM grid with neurons colored based on the value of the feature being plotted.
    """

        # Color scheme
        cmap = mcolors.LinearSegmentedColormap.from_list(
            name='cmap',
            colors=[(0, '#008000'), (0.5, '#FFFF00'), (1, '#FF0000')]
        )

        # Create plot for each feature
        for feature_idx in range(self.train_dat.shape[1]):
            fig, ax = self.plot_map_grid()
            feature_weights = self.weights.iloc[:, feature_idx]
            feature_weights_range = (np.min(feature_weights), np.max(feature_weights))

            # Fill grid with neuron values
            for neuron_idx, row in self.neuron_coordinates.iterrows():
                neuron_value = feature_weights[neuron_idx]
                neuron_fill_color = SOM._map_value_to_color(
                    value=neuron_value,
                    value_range=feature_weights_range,
                    cmap=cmap
                )
                SOM._draw_circle(
                    ax=ax,
                    x=row['x'],
                    y=row['y'],
                    fill_color=neuron_fill_color,
                    edge_color='none'
                )

            SOM._add_colorbar(
                fig=fig,
                value_range=feature_weights_range,
                cmap=cmap
            )
            fig.suptitle(f"Feature {feature_idx}", y=0.95)
            plt.savefig(os.path.join(output_dir, f"feature{feature_idx}_component_plane.png"))
            plt.close(fig)


    def plot_categorical_data(
        self,
        output_dir: str
    ):
        """
        Generate and save categorical data distribution plots across the SOM grid.

        This method visualizes how the SOM, trained using `train_dat`, organizes the data points
        relative to the separate, categorical dataset (`other_dat`). While the SOM is trained only
        on numerical features in `train_dat`, the resulting clustering/organization of the data
        points may aligns with categories in `other_dat`. By examining these plot, we can explore
        whether and how the organization of the SOM reflects patterns in these categorical
        variables.

        Categorical data plots help in assessing whether the SOM's training on numerical data leads
        to meaningful grouping or separation of associated categorical variables. These plots can:
            - Reveal clusters of similar categories across the SOM grid.
            - Provide insights into the relationship between numerical and categorical data.
            - Suggest whether the SOM's organization inherently captures distinctions present in
              the categorical dataset.

        Args:
            - output_dir (str): The directory where the categorical data distribution plots will be
              saved.

        Saves:
            - One plot per categorical feature as a .png file in the specified directory. Each plot displays
              the SOM grid with neurons colored based on the value of the feature being plotted.
    """

        for feature_idx in range(self.other_dat.shape[1]):

            # Initialize the figure and axis for the SOM map grid.
            fig, ax = self.plot_map_grid()

            # Get unique categories and their associated colors.
            categories = np.unique(self.other_dat[:, feature_idx])
            category_colors = SOM._get_distinct_colors(categories)

            # Calculate category proportions for each neuron.
            category_neuron = pd.DataFrame({
                "Neuron": self.observation_mapping,
                "Category": self.other_dat[:, feature_idx]
            })
            counts = category_neuron.groupby(['Neuron', 'Category']).size().unstack(fill_value=0)
            counts = counts.reindex(self.neuron_coordinates.index, fill_value=0)
            category_proportions = counts.div(counts.sum(axis=1), axis=0)

            # Plot pie charts at each neuron's coordinates.
            for idx, (x, y) in self.neuron_coordinates[['x', 'y']].iterrows():

                proportions = category_proportions.loc[idx]
                colors = [category_colors[cat] for cat in proportions.index]

                ax.pie(
                    proportions,
                    colors=colors,
                    center=(x, y),
                    radius=0.5,
                    wedgeprops=dict(edgecolor='black')
                )
            ax.set_aspect('equal')
            ax.set_xlim([-1.5, self.xdim])
            ax.set_ylim([-1, self.ydim - 0.5])

            # Save the plot with a title for the current feature.
            fig.suptitle(f"Feature {feature_idx}", y=0.95)
            output_path = os.path.join(output_dir, f"Feature{feature_idx}_category_proportions.png")
            plt.savefig(output_path)
            plt.close(fig)


    def plot_map_grid(
        self,
        print_neuron_idx: bool = False
    ) -> Tuple[plt.Figure, plt.Axes]:
        """
        Plot a blank SOM grid with optional neuron ids
        """

        fig, ax = plt.subplots(figsize=(self.xdim, self.ydim))

        # Plot each neuron as a circle
        for idx, row in self.neuron_coordinates.iterrows():
            x, y = row['x'], row['y']
            SOM._draw_circle(
                ax=ax,
                x=x,
                y=y,
                fill_color='none',
                edge_color='black'
            )
            # Show neuron id in each neuron
            if print_neuron_idx:
                ax.text(
                    x=x,
                    y=y,
                    s=idx,
                    ha='center',
                    va='center',
                    fontsize=20,
                    color='black'
                )

        ax.set_aspect('equal')
        ax.set_xlim([-1.5, self.xdim])
        ax.set_ylim([-1, self.ydim - 0.5])
        ax.axis('off')
        fig.subplots_adjust(
            left=0.05,
            right=0.95,
            top=0.95,
            bottom=0.05
        )

        return fig, ax

    @staticmethod
    def _draw_circle(
        ax,
        x: float,
        y: float,
        fill_color: Union[Tuple[float, float, float, float], str],
        edge_color: Union[Tuple[float, float, float, float], str],
        radius: float = 0.5
    ):
        """
        Draw a circle on the provided ax object at specified coordinates
        """
        circle = patches.Circle(
                    xy=(x, y),
                    radius=radius,
                    facecolor=fill_color,
                    edgecolor=edge_color
                )
        ax.add_patch(circle)

    @staticmethod
    def _map_value_to_color(
        value: float,
        value_range: Tuple[float, float],
        cmap: mcolors.LinearSegmentedColormap
    ) -> Tuple[float, float, float, float]:

        # Normalize the value to the range [0, 1]
        normalized_value = (value - value_range[0]) / (value_range[1] - value_range[0])

        # Map the normalized value to a color in the gradient
        return cmap(normalized_value)

    @staticmethod
    def _add_colorbar(
        fig: matplotlib.figure.Figure,
        value_range: Tuple[float, float],
        cmap: mcolors.LinearSegmentedColormap
    ):
        """
        Add colorbar to component plane figures
        """
        norm = plt.Normalize(
            vmax=value_range[1],
            vmin=value_range[0]
        )
        sm = plt.cm.ScalarMappable(
            cmap=cmap,
            norm=norm
        )
        sm.set_array([])

        fig.subplots_adjust(bottom=0.1)
        cbar_ax = fig.add_axes([0.125, 0.1, 0.775, 0.03])
        cbar = fig.colorbar(
            mappable=sm,
            cax=cbar_ax,
            orientation='horizontal'
        )
        cbar.set_label(
            label='Color: Neuron Feature Value',
            fontsize=15,
            labelpad=10
        )
        cbar.ax.tick_params(labelsize=14)


    @staticmethod
    def _get_distinct_colors(
        categories: List[str]
    ) -> Dict[str, str]:
        """
        Get a list of n distinct colors using a seaborn color palette.

        Arguements:
            Categories (List[str]): The number of colors required.

        Returns:
            Dict[str, str]: A dictionary mapping each category to a color.
        """
        palette = sns.color_palette('tab20', len(categories))
        palette = palette.as_hex()

        return {categories[i]: palette[i] for i in range(len(categories))}
