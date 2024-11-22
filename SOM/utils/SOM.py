"""
This module define Self-Organize Map (SOM) with asociated attributions, 
mainly training, gettting obseravation, and plotting SOM graph.
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
    """ Building SOM Architecture and its methods.
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
        """Initialize the SOM object with the provided parameters.

        Args:
            train_dat (Union[pd.DataFrame, np.ndarray]): Data to train the SOM on
            other_dat (Union[pd.DataFrame, np.ndarray]): Additional data used to plot on the SOM
            scale_method (str): Method to scale the data. Either 'zscore' or 'minmax'
            x_dim (int): X dimension of the SOM grid
            y_dim (int): Y dimension of the SOM grid
            topology (str): Topology of the SOM grid. Either 'rectangular' or 'hexagonal'
            neighborhood_fnc (str): Neighborhood function to use during training
            epochs (int): Number of epochs to train the SOM for
        """
        # converting Panda dataframe into Numpy array 
        self.train_dat = train_dat.to_numpy() if isinstance(train_dat, pd.DataFrame) else train_dat
        self.other_dat = other_dat.to_numpy() if isinstance(other_dat, pd.DataFrame) else other_dat

        # Scaling and unscaling methods
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
        """ Train the SOM and extract neuron coordinates, weights, and observation mappings

        Args:
            self: SOM object with training data and parameters
        """
        n_samples = self.train_dat_scaled.shape[0]
        n_features = self.train_dat_scaled.shape[1]

        som = MiniSom(
            x=self.xdim,
            y=self.ydim,
            input_len=n_features,
            sigma=1, #related to the how many neighbors are considered
            learning_rate=0.5, # define how much weight is adjusted
            neighborhood_function=self.neighborhood_fnc, # defined neighborhood function
            topology=self.topology,
            activation_distance='euclidean', # method for distance calculation
            random_seed=0
        )

        # principle_component_analysis for weight initialization
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
        """ Z_score calculation: (orignal data - mean) / (standard deviation)
        The mean and standard deviation from original data are used for scaling.
        
        Args:
            self: original training data in numpy array
        
        Returns:
            scaled data (z_score), [mean, standard deviation]: mean and standard deviation of the original data
        """
        means = np.mean(self.train_dat, axis=0)
        stds = np.std(self.train_dat, axis=0)
        return (self.train_dat - means) / stds, [means, stds]


    def _minmax_scale(self):
        """ minmax_score calculation: (orignal data - min) / (max - min)
        It determine minimum and maximum values from the data and cooperate that into scaling.
        
        Args:
            self: original training data in numpy array
        
        Returns:
            scaled data (minmax_score), [min, max]: minimum and maximum values of the original data
        """
        min_vals = np.min(self.train_dat, axis=0)
        max_vals = np.max(self.train_dat, axis=0)
        return (self.train_dat - min_vals) / (max_vals - min_vals), [min_vals, max_vals]


    def _zscore_unscale(self, scaled_data):
        """ Reverse scaled to unscaled data z_score

        Args:
            scaled_data (numpy array): scaled data

        Returns:
            unscaled data (z_score): original data
        """
        means, stds = self._scaling_factors
        return scaled_data * stds + means


    def _minmax_unscale(self, scaled_data):
        """ Reverse scaled to unscaled data minmax_score

        Args:
            scaled_data (numpy array): scaled data

        Returns:
            unscaled data (minmax_score): original data
        """
        min_vals, max_vals = self._scaling_factors
        return scaled_data * (max_vals - min_vals) + min_vals


    def _get_observation_neuron_mappings(self):
        """ Get the neuron id of the winning neuron for each observation in the training data

        Returns:
            List[int]: List of neuron ids for each observation
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
        """ Get the x-y coordinates of each neuron in the SOM grid

        Returns:
            pd.DataFrame: DataFrame with columns 'x' and 'y' containing the coordinates
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
        """ Get the weights of each neuron in the SOM grid

        Returns:
            pd.DataFrame: DataFrame with columns for each feature and rows for each neuron
        """
        weights_scaled = [
            self.map.get_weights()[i, j, :] for i in range(self.xdim) for j in range(self.ydim)
        ]
        return pd.DataFrame(weights_scaled)


    def plot_component_planes(
        self,
        output_dir: str,
        save: bool = True
    ):
        """ Plot the component planes of the SOM

        Args:
            output_dir (str): Directory to save the component plane plots
            save (bool): Whether to save the plots to disk, default is True
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
            if save:
                plt.savefig(os.path.join(output_dir, f"feature{feature_idx}_component_plane.png"))
            plt.close(fig)


    def plot_categorical_data(
        self,
        output_dir: str,
        save: bool = True
    ):
        """ Plot the categorical data on the SOM grid

        Args:
            output_dir (str): Directory to save the plots
            save (bool): Whether to save the plots to disk, default is True
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
            if save:
                plt.savefig(output_path)
            plt.close(fig)


    def plot_map_grid(
        self,
        print_neuron_idx: bool = False
    ) -> Tuple[plt.Figure, plt.Axes]:
        """ Plot the SOM grid with each neuron as a circle

        Args:
            print_neuron_idx (bool, optional): Whether to print the neuron index in each neuron. Defaults to False.

        Returns:
            Tuple[plt.Figure, plt.Axes]: Figure and axis objects for the plot
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
        """ Draw a circle on the plot

        Args:
            ax: Axis object to draw the circle on
            x (float): X coordinate of the circle
            y (float): Y coordinate of the circle
            fill_color (Union[Tuple[float, float, float, float], str]): Fill color of the circle
            edge_color (Union[Tuple[float, float, float, float], str]): Edge color of the circle
            radius (float, optional): Radius of the circle. Defaults to 0.5.
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
        """ Map a value to a color in a given colormap

        Args:
            value (float): value to map to a color
            value_range (Tuple[float, float]): Range of values to normalize the value to
            cmap (mcolors.LinearSegmentedColormap): Colormap to map the value to

        Returns:
            Tuple[float, float, float, float]: RGBA color tuple
        """

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
        """ Add a colorbar to the figure

        Args:
            fig (matplotlib.figure.Figure): Figure to add the colorbar to
            value_range (Tuple[float, float]): Range of values to normalize the colorbar to
            cmap (mcolors.LinearSegmentedColormap): Colormap to use for the colorbar
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