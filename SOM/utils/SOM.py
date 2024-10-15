import os
from typing import Tuple, Union

import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import matplotlib.figure
from minisom import MiniSom
import numpy as np
import pandas as pd


class SOM():
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
        """_
        Train the SOM and extract neuron coordinates, weights, and observation mappings
        """
        n_samples = self.train_dat_scaled.shape[0]
        n_features = self.train_dat_scaled.shape[1]

        som = MiniSom(
            x=self.xdim,
            y=self.ydim,
            input_len=n_features,
            sigma=1,
            learning_rate=0.5,
            neighborhood_function=self.neighborhood_fnc,
            topology=self.topology,
            activation_distance='euclidean',
            random_seed=0
        )
        som.pca_weights_init(self.train_dat_scaled)
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
        means = np.mean(self.train_dat, axis=0)
        stds = np.std(self.train_dat, axis=0)
        return (self.train_dat - means) / stds, [means, stds]

    def _minmax_scale(self):
        min_vals = np.min(self.train_dat, axis=0)
        max_vals = np.max(self.train_dat, axis=0)
        return (self.train_dat - min_vals) / (max_vals - min_vals), [min_vals, max_vals]

    def _zscore_unscale(self, scaled_data):
        means, stds = self._scaling_factors
        return scaled_data * stds + means

    def _minmax_unscale(self, scaled_data):
        min_vals, max_vals = self._scaling_factors
        return scaled_data * (max_vals - min_vals) + min_vals

    def _get_observation_neuron_mappings(self):
        """
        Get the neuron id that each observation in the training data is assigned to
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
        Get the x-y coordinates of the SOM grid's neurons
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
        Get the weights of each neuron after training on scaled data.
        """
        # weights_scaled = []  # "scaled" since trained on scaled data
        # for i in range(self.x_dim * self.ydim):
        #     neuron_idx = np.unravel_index(i, (self.x_dim, self.ydim))
        #     row, col = neuron_idx[0], neuron_idx[1]
        #     weight_scaled = som.get_weights()[row, col, :]
        #     weights_scaled.append(weight_scaled)
        # weights_scaled = pd.DataFrame(weights_scaled)

        weights_scaled = [
            self.map.get_weights()[i, j, :] for i in range(self.xdim) for j in range(self.ydim)
        ]
        return pd.DataFrame(weights_scaled)

    def plot_component_planes(
        self,
        output_dir: str
    ):
        """
        Plot component planes for each feature in the training data and save the figures.
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
