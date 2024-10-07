from typing import Union
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
        self.train_dat = train_dat if isinstance(train_dat, np.ndarray) else train_dat.to_numpy()
        self.other_dat = other_dat if isinstance(other_dat, np.ndarray) else other_dat.to_numpy()
        self.scale = getattr(self, f"{scale_method}_scale", None)
        self.x_dim = x_dim
        self.ydim = y_dim
        self.topology = topology
        self.neighborhood_fnc = neighborhood_fnc
        self.epochs = epochs

    def train_map(self):
        train_dat_scaled = self.scale()
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
        self.map = som

    def zscore_scale(self):
        mean = np.mean(self.train_dat, axis=0)
        std = np.std(self.train_dat, axis=0)
        return (self.train_dat - mean) / std

    def minmax_scale(self):
        min_val = np.min(self.train_dat, axis=0)
        max_val = np.max(self.train_dat, axis=0)
        return (self.train_dat - min_val) / (max_val - min_val)

    def plot_component_planes(self):
        return
    
    def plot_categories(self):
        return
    
    def plot_map_grid(self):
        return
    