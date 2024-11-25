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

import math
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
    A class for building, training, and visualizing Self-Organizing Maps (SOMs).

    This class provides tools to train a Self-Organizing Map on high-dimensional data, including
    scaling and unscaling data, calculating various fit metrics (e.g., topographic error, percent
    variance explained), and visualizing the SOM grid after fitting.

    Attributes (set before training):
    ---------------------------------
    train_dat (np.ndarray): 
        The raw (unscaled) training data used to fit the SOM.
    other_dat (np.ndarray): 
        Additional, non-numeric data associated with each data point in `train_dat`, but 
        not used to train the SOM.
    train_dat_features (List[str]): 
        A list of feature names associated with `train_dat`.
    other_dat_features (List[str]): 
        A list of feature names associated with `other_dat`.
    xdim (int): 
        The number of neurons in the x-dimension of the SOM grid.
    ydim (int): 
        The number of neurons in the y-dimension of the SOM grid.
    topology (str): 
        The topology of the SOM grid ('rectangular' or 'hexagonal', associated with 4 and 6
        neighbors, respectively).
    neighborhood_fnc (str): 
        The neighborhood function used during training ('gaussian' or 'bubble').
    epochs (int): 
        The number of training epochs (i.e., how many times the dataset is presented during the
        training process).
    _scale (Callable): 
        Method for scaling `train_dat` (z-score or min/max).
    _unscale (Callable): 
        Method for unscaling scaled data (z-score or min/max).

    Attributes (set after training):
    --------------------------------
    map (MiniSom): 
        The trained MiniSom instance, created after calling `train_map`.
    observation_mapping (np.ndarray): 
        The mapping of data points to their closest neurons (i.e. BMU, best matching unit).
    neuron_coordinates (pd.DataFrame): 
        The coordinates of neurons in the SOM grid.
    weights (pd.DataFrame): 
        The unscaled weight vectors of each neuron after training.
    weights_scaled (pd.DataFrame): 
        The scaled weight vectors of each neuron after training.

    Methods:
    --------
    __init__:
        Initializes the SOM instance with training data, configuration parameters, and scaling 
        methods.

    _validate_inputs:
        Validates the input parameters for the SOM class, ensuring data types, ranges, and formats 
        are correct.

    train_map:
        Trains the SOM using scaled training data. Includes initializing neuron weights and mapping
        data points to neurons.

    calculate_percent_variance_explained:
        Calculates the Percent Variance Explained (PVE) of the fit SOM, indicating how well the
        SOM represents the variance in the data.

    calculate_topographic_error:
        Calculates the topographic error of the fit SOM, a measure of how well the SOM preserves
        the topology of the data.

    plot_component_planes:
        Generates and saves component plane plots for each feature in the training data.

    plot_categorical_data:
        Generates and saves the distribution of categorical data (from `other_dat`) across the 
        SOM grid.

    plot_map_grid:
        Plots a blank SOM grid, optionally labeling neurons with their indices.

    _zscore_scale:
        Scales the training data (`train_dat`) using z-score normalization.

    _minmax_scale:
        Scales the training data (`train_dat`) using min-max normalization.

    _zscore_unscale:
        Reverses z-score normalization to unscale data back to its original range.

    _minmax_unscale:
        Reverses min-max scaling to unscale data back to its original range.

    _get_observation_neuron_mappings:
        Maps each data point in the training set (`train_dat`) to its Best Matching Unit (BMU).

    _get_neuron_coordinates:
        Retrieves the x and y coordinates of all neurons in the SOM grid.

    _get_weights:
        Retrieves the weight vectors of all neurons in the SOM grid after training.

    _draw_circle:
        Draws a circle at specified coordinates in a matplotlib plot, representing a neuron.

    _map_value_to_color:
        Maps numeric feature values to colors using a specified colormap.

    _add_colorbar:
        Adds a colorbar to a matplotlib plot to represent feature value ranges.

    _get_distinct_colors:
        Generates distinct colors for categorical data using a seaborn palette.

    _check_output_path:
        Ensures the specified directory exists, creating it if necessary.

    _are_nodes_adjacent:
        Determines if two neurons in the SOM grid are adjacent based on their coordinates, used
        when calculating topographic error.
    """

    def __init__(
        self,
        train_dat: pd.DataFrame,
        other_dat: pd.DataFrame,
        scale_method: str,
        x_dim: int,
        y_dim: int,
        topology: str,
        neighborhood_fnc: str,
        epochs: int
    ):

        """
        Initialize the Self-Organizing Map (SOM) class with user-defined parameters.

        This constructor sets up the SOM by validating inputs, scaling the training data, and
        preparing the attributes needed for training and analysis.

        Args:
            train_dat (pd.DataFrame): A DataFrame containing the numerical training data for the 
                SOM.
            other_dat (pd.DataFrame): A DataFrame with non-numeric data associated with each row in
                `train_dat`, not used for training the SOM but for analysis.
            scale_method (str): Method for scaling the data; must be either 'zscore' or 'minmax'.
            x_dim (int): The number of neurons in the x-dimension of the SOM grid. Must be a 
                positive integer.
            y_dim (int): The number of neurons in the y-dimension of the SOM grid. Must be a
                positive integer.
            topology (str): Topology of the SOM grid; must be 'rectangular' or 'hexagonal'.
            neighborhood_fnc (str): Neighborhood function to use during training; must be
                'gaussian' or 'bubble'.
            epochs (int): Number of epochs (iterations over the dataset) for training the SOM. Must
                be a positive integer.
        """

        # Input validation
        self._validate_inputs(
            train_dat=train_dat,
            other_dat=other_dat,
            scale_method=scale_method,
            x_dim=x_dim,
            y_dim=y_dim,
            topology=topology,
            neighborhood_fnc=neighborhood_fnc,
            epochs=epochs
        )

        # Save feature names and convert input data to numpy array
        self.train_dat_features = train_dat.columns.tolist()
        self.other_dat_features = other_dat.columns.tolist()
        self.train_dat = train_dat.to_numpy()
        self.other_dat = other_dat.to_numpy()

        # Scaling and attributes
        self._scale = getattr(self, f"_{scale_method}_scale", None)
        self._unscale = getattr(self, f"_{scale_method}_unscale", None)
        self.train_dat_scaled, self._scaling_factors = self._scale()
        self.xdim = x_dim
        self.ydim = y_dim
        self.topology = topology
        self.neighborhood_fnc = neighborhood_fnc
        self.epochs = epochs

        # Post-training attributes
        self.map = None
        self.observation_mapping = None
        self.neuron_coordinates = None
        self.weights = None
        self.weights_scaled = None


    @staticmethod
    def _validate_inputs(
        train_dat: pd.DataFrame,
        other_dat: pd.DataFrame,
        scale_method: str,
        x_dim: int,
        y_dim: int,
        topology: str,
        neighborhood_fnc: str,
        epochs: int):

        """
        Validate the input parameters for initializing the SOM class.

        Ensures that all required inputs meet the expected types, ranges, and formats.

        Args:
            train_dat (pd.DataFrame): A DataFrame containing the numerical training data for the 
                SOM.
            other_dat (pd.DataFrame): A DataFrame with non-numeric data associated with each row in
                `train_dat`, not used for training the SOM but for analysis.
            scale_method (str): Method for scaling the data; must be either 'zscore' or 'minmax'.
            x_dim (int): The number of neurons in the x-dimension of the SOM grid. Must be a 
                positive integer.
            y_dim (int): The number of neurons in the y-dimension of the SOM grid. Must be a
                positive integer.
            topology (str): Topology of the SOM grid; must be 'rectangular' or 'hexagonal'.
            neighborhood_fnc (str): Neighborhood function to use during training; must be
                'gaussian' or 'bubble'.
            epochs (int): Number of epochs (iterations over the dataset) for training the SOM. Must
                be a positive integer.
            
        Raises:
            TypeError: If `train_dat` or `other_dat` is not a pandas DataFrame, or if `train_dat`
                contains non-numeric data.
            ValueError: If `scale_method`, `topology`, or `neighborhood_fnc` contains invalid
                values, or if `x_dim`, `y_dim`, or `epochs` are not positive integers.
        """

        if not isinstance(train_dat, pd.DataFrame):
            raise TypeError("train_dat must be a pandas DataFrame.")

        is_all_numeric = train_dat.dtypes.apply(lambda dtype: np.issubdtype(dtype, np.number)).all()
        if not is_all_numeric:
            raise TypeError("train_dat must contain only numeric values.")

        if not isinstance(other_dat, pd.DataFrame):
            raise TypeError("other_dat must be a pandas DataFrame.")

        if scale_method not in ["zscore", "minmax"]:
            raise ValueError("scale_method must be 'zscore' or 'minmax'.")

        if topology not in ["rectangular", "hexagonal"]:
            raise ValueError("topology must be 'rectangular' or 'hexagonal'.")

        if neighborhood_fnc not in ["gaussian", "bubble"]:
            raise ValueError("neighborhood_fnc must be 'gaussian' or 'bubble'.")

        if not isinstance(x_dim, int) or x_dim <= 0:
            raise ValueError("x_dim must be a positive integer.")

        if not isinstance(y_dim, int) or y_dim <= 0:
            raise ValueError("y_dim must be a positive integer.")

        if not isinstance(epochs, int) or epochs <= 0:
            raise ValueError("epochs must be a positive integer.")


    def train_map(self):

        """
        Train the Self-Organizing Map on the scaled training data.

        This method initializes the SOM using PCA-based weight initialization and trains it using
        the specified parameters. It also calculates neuron weights and maps observations to their
        corresponding neurons.

        Sets the following attributes:
            - map: The trained MiniSom instance.
            - observation_mapping: Mapping of observations to neurons.
            - neuron_coordinates: Coordinates of neurons in the SOM grid.
            - weights_scaled: Scaled weights of each neuron.
            - weights: Unscaled weights of each neuron.
        """

        n_samples = self.train_dat_scaled.shape[0]
        n_features = self.train_dat_scaled.shape[1]

        if n_features <= 1:
            raise ValueError("Training data must have at least two features.")

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
        Get the neuron id that each observation in `train_dat` is assigned to.

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
        Retrieve the weight vector of each neuron after training.

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
        in `train_dat` are distributed across the SOM grid. Each neuron on the grid is assigned
        a color based on the value of the feature it represents, allowing you to see patterns and
        relationships in the data.

        Component planes help in understanding how individual features contribute to the
        organization of the SOM. By examining these plots, you can:
            - Identify clusters and trends in specific features.
            - Compare how different features vary across the grid, potentially highlighting
                relationships/correlations between features.

        Args:
            output_dir (str): The directory where the component plane plots will be saved. 

        Saves:
            One plot per feature as a .png file in the specified directory. Each plot displays
            the SOM grid with neurons colored based on the value of the feature being plotted.
        """

        # Check that SOM has been trained
        if self.map is None:
            raise RuntimeError("SOM has not been trained. Call `train_map` before plotting.")

        # Color scheme
        cmap = mcolors.LinearSegmentedColormap.from_list(
            name='cmap',
            colors=[(0, '#008000'), (0.5, '#FFFF00'), (1, '#FF0000')]
        )

        # Create one plot for each feature
        for feature_idx in range(self.train_dat.shape[1]):
            fig, ax = self.plot_map_grid()
            feature_weights = self.weights.iloc[:, feature_idx]
            feature_weights_range = (np.min(feature_weights), np.max(feature_weights))
            feature_name = self.train_dat_features[feature_idx]

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
                    edge_color='none',
                    text=str(round(neuron_value, 2))
                )

            SOM._add_colorbar(
                fig=fig,
                value_range=feature_weights_range,
                cmap=cmap,
                label=feature_name
            )
            fig.suptitle(f"{feature_name}", y=1)

            fig.subplots_adjust(
                top=0.9,
                bottom=0.15,
                left=0.1,
                right=0.9
            )

            # Set output directory and create if it doesn't exist
            output_path = os.path.join(output_dir, f"{feature_name}_component_plane.png")
            self._check_output_path(
                path=output_path
            )

            plt.savefig(
                fname=output_path,
                format='png',
                dpi=300,
                bbox_inches='tight')
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
            output_dir (str): The directory where the categorical data distribution plots will be
                saved.

        Saves:
            One plot per categorical feature as a .png file in the specified directory. Each plot
            displays the SOM grid with neurons colored based on the value of the feature being
            plotted.
        """

        # Check that SOM has been trained
        if self.map is None:
            raise RuntimeError("SOM has not been trained. Call `train_map` before plotting.")

        # Make one plot per categorical feature
        for feature_idx in range(self.other_dat.shape[1]):

            feature_name = self.other_dat_features[feature_idx]

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

            # Add legend
            legend_patches = [
                patches.Patch(color=color, label=cat) for cat, color in category_colors.items()
            ]
            ax.legend(
                handles=legend_patches,
                loc='upper center',
                bbox_to_anchor=(0.5, 0.1),
                ncol=3,
                fontsize=9,
                frameon=False
            )

            # Add title
            fig.suptitle(f"Proportion of {feature_name} per neuron", y=1)

            fig.subplots_adjust(
                top=0.9,
                bottom=0.2,
                left=0.1,
                right=0.9
            )

            # Set output directory and create if it doesn't exist
            output_path = os.path.join(output_dir, f"{feature_name}_category_proportions.png")
            self._check_output_path(
                path=output_path
            )

            plt.savefig(
                fname=output_path,
                format='png',
                dpi=300,
                bbox_inches='tight'
            )
            plt.close(fig)


    def plot_map_grid(
        self,
        print_neuron_idx: bool = False
    ) -> Tuple[plt.Figure, plt.Axes]:

        """ 
        Plot a blank SOM grid where each neuron is represented by a circle.

        This method creates a visual representation of the SOM grid without any feature or category
        values. Optionally, you can display the index of each neuron within its circle.

        Args:
            print_neuron_idx (bool, optional): Whether to display the neuron indices on the plot. 
                Defaults to False.

        Returns:
            Tuple[plt.Figure, plt.Axes]: The matplotlib figure and axis objects for the plot.
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
        ax: plt.axes,
        x: float,
        y: float,
        fill_color: Union[Tuple[float, float, float, float], str],
        edge_color: Union[Tuple[float, float, float, float], str],
        text: str = None,
        radius: float = 0.5
    ):

        """ 
        Draw a circle at the specified coordinates on a matplotlib axis.

        Args:
            ax (matplotlib.axes.Axes): The axis object where the circle will be drawn.
            x (float): The x-coordinate of the circle's center.
            y (float): The y-coordinate of the circle's center.
            fill_color (Union[Tuple[float, float, float, float], str]): The color to fill the
                circle with.
            edge_color (Union[Tuple[float, float, float, float], str]): The color of the circle's
                edge.
            text (str, optional): Text to annotate neuron.
            radius (float, optional): The radius of the circle. Defaults to 0.5.
        """

        circle = patches.Circle(
                    xy=(x, y),
                    radius=radius,
                    facecolor=fill_color,
                    edgecolor=edge_color
                )
        ax.add_patch(circle)

        if text:
            ax.text(
                x=x,
                y=y,
                s=text,
                ha='center',
                va='center',
                fontsize=8,
                color='black'
                )


    @staticmethod
    def _map_value_to_color(
        value: float,
        value_range: Tuple[float, float],
        cmap: mcolors.LinearSegmentedColormap
    ) -> Tuple[float, float, float, float]:

        """ 
        Map a numeric value to a color using a specified colormap.

        Args:
            value (float): The numeric value to map.
            value_range (Tuple[float, float]): The range of possible values (min, max).
            cmap (mcolors.LinearSegmentedColormap): The colormap to use for mapping.

        Returns:
            Tuple[float, float, float, float]: The RGBA color corresponding to the input value.
        """

        # Normalize the value to the range [0, 1]
        normalized_value = (value - value_range[0]) / (value_range[1] - value_range[0])

        # Map the normalized value to a color in the gradient
        return cmap(normalized_value)


    @staticmethod
    def _add_colorbar(
        fig: matplotlib.figure.Figure,
        value_range: Tuple[float, float],
        cmap: mcolors.LinearSegmentedColormap,
        label: str
    ):

        """ 
        Add a horizontal colorbar to a matplotlib figure.

        Args:
            fig (matplotlib.figure.Figure): The figure object where the colorbar will be added.
            value_range (Tuple[float, float]): The range of values (min, max) represented by the
                colorbar.
            cmap (mcolors.LinearSegmentedColormap): The colormap to use for the colorbar.
            label (str): Name of feature
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
        cbar_ax = fig.add_axes([0.2, 0.2, 0.6, 0.02])
        cbar = fig.colorbar(
            mappable=sm,
            cax=cbar_ax,
            orientation='horizontal'
        )
        cbar.set_label(
            label=f'Color: Neuron {label} Value',
            fontsize=9,
            labelpad=10
        )
        cbar.ax.tick_params(labelsize=9)


    @staticmethod
    def _get_distinct_colors(
        categories: List[str]
    ) -> Dict[str, str]:

        """
        Generate a dictionary of distinct colors for a list of categories.

        This method uses the seaborn color palette to assign a unique color to each category.

        Args:
            categories (List[str]): A list of unique category labels.

        Returns:
            Dict[str, str]: A dictionary mapping each category to a distinct hex color string.
        """

        palette = sns.color_palette('tab20', len(categories))
        palette = palette.as_hex()

        return {categories[i]: palette[i] for i in range(len(categories))}


    @staticmethod
    def _check_output_path(
        path: str
    ):
        
        """
        Ensure the specified output path exists. If the directory does not exist, it is created.

        Args:
            path (str): The directory path to check or create.

        Raises:
            OSError: If the directory cannot be created due to permissions or other errors.
        """

        if not os.path.exists(path):
            os.makedirs(path)


    def calculate_topographic_error(self) -> float:

        """
        Calculate the topographic error of the trained SOM.

        Topographic error measures the proportion of data points for which the two closest neurons
        (BMUs) are not adjacent in the SOM grid. A lower topographic error indicates better 
        topology preservation.

        Returns:
            float: The topographic error, ranging from 0 (perfect preservation) to 1 (poor 
                preservation).
        """

        if self.map is None:
            raise RuntimeError(
                "SOM has not been trained. Call `train_map` before calculating topographic error."
            )

        topographic_error_count = 0

        for idx in np.arange(len(self.train_dat_scaled)):
            # Find the BMU (Best Matching Unit) for the current data point
            bmu1 = self.observation_mapping[idx]

            # Get the weights for all nodes
            all_weights = self.weights_scaled.to_numpy()

            # Compute distances from current data point to all nodes
            distances = np.linalg.norm(all_weights - self.train_dat_scaled[idx, :], axis=1)

            # Identify the second-best matching unit
            bmu2 = np.argsort(distances)[1]

            # Check if the two BMUs are adjacent in the SOM grid
            if not self._are_nodes_adjacent(bmu1, bmu2):
                topographic_error_count += 1

        # Calculate the topographic error as the proportion of non-adjacent BMUs
        return topographic_error_count / len(self.train_dat_scaled)


    def calculate_percent_variance_explained(self) -> float:

        """
        Calculate the Percent Variance Explained (PVE) for the trained SOM.

        PVE is the proportion of the total variance in the data explained by the SOM.
        It is computed as: PVE = ((TSS - WCSS) / TSS) * 100, where:
            - TSS = total sum of squares
            - WCSS = within cluster (i.e., within neuron) sum of squares

        Returns:
            float: Percent Variance Explained (PVE)
        """

        if self.map is None:
            raise RuntimeError(
                "SOM has not been trained. Call `train_map` before calculating PVE."
            )

        # Compute the mean of the dataset
        global_mean = self.train_dat_scaled.mean(axis=0)

        # Compute Total Sum of Squares (TSS)
        tss = np.sum(np.linalg.norm(self.train_dat_scaled - global_mean, axis=1) ** 2)

        # Compute Within-Cluster Sum of Squares (WCSS)
        wcss = 0
        for idx in np.arange(len(self.train_dat_scaled)):
            bmu = self.observation_mapping[idx]
            bmu_weight = self.weights_scaled.iloc[bmu, :].to_numpy()
            wcss += np.linalg.norm(self.train_dat_scaled[idx, :] - bmu_weight) ** 2

        # Calculate PVE
        pve = ((tss - wcss) / tss) * 100

        return pve


    def _are_nodes_adjacent(
        self,
        bmu1: int,
        bmu2: int
    ) -> bool:

        """
        Check if two neurons in the SOM grid are adjacent.

        Args:
            bmu1 (int): ID of the first neuron.
            bmu2 (int): ID of the second neuron.

        Returns:
            bool: True if the neurons are neighbors, False otherwise.
        """

        coord1 = self.neuron_coordinates.iloc[bmu1, :].to_list()
        coord2 = self.neuron_coordinates.iloc[bmu2, :].to_list()

        euclidean_distance = np.sqrt((coord1[0] - coord2[0]) ** 2 + (coord1[1] - coord2[1]) ** 2)

        return math.isclose(
            a=euclidean_distance,
            b=1
        )
