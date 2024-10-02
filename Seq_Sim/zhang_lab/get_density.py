from scipy.stats import gaussian_kde
import numpy as np


def get_density(x, y, bandwidth="scott"):
    # Create a Gaussian KDE object
    kde = gaussian_kde(np.vstack([x, y]), bw_method=bandwidth)

    # Create a grid for evaluation
    x_grid = np.linspace(np.min(x), np.max(x), 100)
    y_grid = np.linspace(np.min(y), np.max(y), 100)
    x_mesh, y_mesh = np.meshgrid(x_grid, y_grid)

    # Evaluate the density on the grid
    z = kde(np.vstack([x_mesh.ravel(), y_mesh.ravel()])).reshape(x_mesh.shape)

    return z
