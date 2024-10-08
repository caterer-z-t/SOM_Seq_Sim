from typing import Union
from minisom import MiniSom
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches


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
        self.train_dat = train_dat if isinstance(train_dat, np.ndarray) else train_dat.to_numpy()
        self.other_dat = other_dat if isinstance(other_dat, np.ndarray) else other_dat.to_numpy()
        self._scale = getattr(self, f"{scale_method}_scale", None)
        self._unscale = getattr(self, f"{scale_method}_unscale", None)
        self.train_dat_scaled, self._scaling_factors = self._scale()
        self.x_dim = x_dim
        self.ydim = y_dim
        self.topology = topology
        self.neighborhood_fnc = neighborhood_fnc
        self.epochs = epochs
        self.map = None
        self.neuron_idx = None
        self.neuron_coordinates = None
        self.weights = None
        self.weights_scaled = None

    def train_map(self):
        train_dat_scaled = self.train_dat_scaled
        n_samples = train_dat_scaled.shape[0]
        n_features = train_dat_scaled.shape[1]

        som = MiniSom(
            x=self.x_dim,
            y=self.ydim,
            input_len=n_features,
            sigma=1,
            learning_rate=0.5,
            neighborhood_function=self.neighborhood_fnc,
            topology=self.topology,
            activation_distance='euclidean',
            random_seed=0
        )
        som.train(
            data=train_dat_scaled,
            num_iteration=self.epochs * n_samples,
            random_order=True,
            verbose=False
        )
        # Map each observation to neuron
        winner_coordinates = np.array([som.winner(x) for x in train_dat_scaled]).T
        neuron_idx = np.ravel_multi_index(winner_coordinates, (self.x_dim, self.ydim))

        # Get coordiantes of each neuron
        xx, yy = som.get_euclidean_coordinates()
        weights = som.get_weights()

        x_coords = [xx[(i, j)] for i in range(weights.shape[0]) for j in range(weights.shape[1])]
        y_coords = [yy[(i, j)] for i in range(weights.shape[0]) for j in range(weights.shape[1])]
        neuron_coordiantes = pd.DataFrame({'x': x_coords, 'y': y_coords})

        # Get weights of each neuron
        weights_scaled = []  # "scaled" since trained on scaled data
        for i in range(self.x_dim * self.ydim):
            neuron_idx = np.unravel_index(i, (self.x_dim, self.ydim))
            row, col = neuron_idx[0], neuron_idx[1]
            weight_scaled = som.get_weights()[row, col, :]
            weights_scaled.append(weight_scaled)
        weights_scaled = pd.DataFrame(weights_scaled)

        # Set attributes
        self.map = som
        self.train_dat_scaled = train_dat_scaled
        self.neuron_idx = neuron_idx
        self.neuron_coordinates = neuron_coordiantes
        self.weights_scaled = weights_scaled
        self.weights = self._unscale(weights_scaled)

    def zscore_scale(self):
        means = np.mean(self.train_dat, axis=0)
        stds = np.std(self.train_dat, axis=0)
        return (self.train_dat - means) / stds, [means, stds]

    def minmax_scale(self):
        min_vals = np.min(self.train_dat, axis=0)
        max_vals = np.max(self.train_dat, axis=0)
        return (self.train_dat - min_vals) / (max_vals - min_vals), [min_vals, max_vals]

    def zscore_unscale(self, scaled_data):
        means, stds = self._scaling_factors
        return scaled_data * stds + means

    def minmax_unscale(self, scaled_data):
        min_vals, max_vals = self._scaling_factors
        return scaled_data * (max_vals - min_vals) + min_vals

    def plot_component_planes(self):
        fig, ax = self.plot_map_grid()

        for 

        return

    def plot_categories(self):
        return

    def plot_map_grid(self, print_neuron_idx=False):

        BASE_SIZE_PER_NEURON = 1

        # Calculate fig size based on SOM dimensions
        fig_width = BASE_SIZE_PER_NEURON * self.x_dim
        fig_height = BASE_SIZE_PER_NEURON * self.ydim

        # Calculate radius of circles (neurons) for plotting
        circle_radius = 0.5

        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

        # Plot each neuron as a circle
        for idx, row in self.neuron_coordinates.iterrows():
            x, y = row['x'], row['y'] * np.sqrt(3) / 2

            circle = patches.Circle(
                xy=(x, y),
                radius=circle_radius,
                edgecolor='black',
                facecolor='none'
            )
            ax.add_patch(circle)

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
        ax.set_xlim([-1.5, self.x_dim])
        ax.set_ylim([-1, self.ydim - 0.5])
        ax.axis('off')
        ax.margins(0)
        fig.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.05)

        return fig, ax
