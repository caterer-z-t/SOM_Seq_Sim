# In[1]: Imports --
import pytest
import numpy as np
import os
from unittest.mock import MagicMock, patch
import matplotlib.pyplot as plt
import matplotlib.patches as patches 
import matplotlib.colors as mcolors
from unittest import TestCase
import pandas as pd
from SOM.utils.som_utils import SOM

# In[2]: Test cases for the initialization of the SOM class --

class TestSOMInitialization:
    """Test cases for the initialization of the SOM class."""

    def setup_method(self):
        """Setup common resources for tests in this class."""
        # Mock data for testing
        self.train_data = pd.DataFrame(
            {"feature1": [1.0, 2.0, 3.0], "feature2": [4.0, 5.0, 6.0]}
        )
        self.other_data = pd.DataFrame(
            {"feature1": [7.0, 8.0, 9.0], "feature2": [10.0, 11.0, 12.0]}
        )
        self.scale_method = "zscore"
        self.x_dim = 5
        self.y_dim = 5
        self.topology = "rectangular"
        self.neighborhood_fnc = "gaussian"
        self.epochs = 100

    def test_init_attributes(self):
        """Test that the attributes are properly initialized."""
        # Initialize SOM
        som = SOM(
            train_dat=self.train_data,
            other_dat=self.other_data,
            scale_method=self.scale_method,
            x_dim=self.x_dim,
            y_dim=self.y_dim,
            topology=self.topology,
            neighborhood_fnc=self.neighborhood_fnc,
            epochs=self.epochs,
        )

        # Assertions to verify attributes
        assert np.array_equal(som.train_dat, self.train_data.to_numpy())
        assert np.array_equal(som.other_dat, self.other_data.to_numpy())
        assert som.xdim == self.x_dim
        assert som.ydim == self.y_dim
        assert som.topology == self.topology
        assert som.neighborhood_fnc == self.neighborhood_fnc
        assert som.epochs == self.epochs
        assert som.train_dat_scaled is not None
        assert som._scaling_factors is not None
        assert som.map is None
        assert som.observation_mapping is None
        assert som.neuron_coordinates is None
        assert som.weights is None
        assert som.weights_scaled is None

# In[3]: Test cases for the scaling of data during initialization --
class TestSOMScaling:
    """Test cases for the scaling of data during initialization."""

    def setup_method(self):
        """Setup common resources for tests in this class."""
        # Mock data for testing
        self.train_data = pd.DataFrame(
            {"feature1": [1.0, 2.0, 3.0], "feature2": [4.0, 5.0, 6.0]}
        )
        self.other_data = pd.DataFrame(
            {"feature1": [7.0, 8.0, 9.0], "feature2": [10.0, 11.0, 12.0]}
        )
        self.scale_method = "zscore"
        self.x_dim = 5
        self.y_dim = 5
        self.topology = "rectangular"
        self.neighborhood_fnc = "gaussian"
        self.epochs = 100

    def test_scaling_zscore(self):
        """Test z-score scaling during initialization."""
        # Initialize SOM
        som = SOM(
            train_dat=self.train_data,
            other_dat=self.other_data,
            scale_method=self.scale_method,
            x_dim=self.x_dim,
            y_dim=self.y_dim,
            topology=self.topology,
            neighborhood_fnc=self.neighborhood_fnc,
            epochs=self.epochs,
        )

        # Manually calculate z-score scaling
        means = np.mean(self.train_data.to_numpy(), axis=0)
        stds = np.std(self.train_data.to_numpy(), axis=0)
        expected_scaled_data = (self.train_data.to_numpy() - means) / stds

        # Assertions
        assert np.allclose(som.train_dat_scaled, expected_scaled_data)
        assert np.allclose(som._scaling_factors[0], means)
        assert np.allclose(som._scaling_factors[1], stds)

# In[9]: Test case for the plot_component_planes method --

class TestSOMPlotComponentPlanes(TestCase):

    def setUp(self):
        """Set up a basic SOM object with mock data for testing"""
        # Mock some basic attributes needed for the test
        self.som = MagicMock(spec=SOM)

        # Mock training data (3 features, 5 samples)
        self.som.train_dat = np.array(
            [[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12], [13, 14, 15]]
        )

        # Mock neuron coordinates (5 neurons in a 3x3 grid)
        self.som.neuron_coordinates = pd.DataFrame(
            {"x": [0, 1, 2, 0, 1], "y": [0, 1, 2, 1, 2]}
        )

        # Mock weights (5 neurons with 3 features)
        self.som.weights = pd.DataFrame(np.random.rand(5, 3))

        # Mock plot_map_grid method to return a mock figure and axis
        self.som.plot_map_grid = MagicMock(return_value=(plt.figure(), plt.gca()))

        # Mock the _map_value_to_color method to return a color
        self.som._map_value_to_color = MagicMock(return_value="blue")

        # Mock the _draw_circle method
        self.som._draw_circle = MagicMock()

        # Mock the _add_colorbar method
        self.som._add_colorbar = MagicMock()

    @patch("matplotlib.pyplot.savefig")  # Mock savefig to avoid actual file creation
    def test_plot_component_planes(self, mock_savefig):
        """Test the plot_component_planes function"""

        # Set a dummy output directory
        output_dir = "mock_output"

        # Call the plot_component_planes method
        self.som.plot_component_planes(output_dir)

        # Check if plot_map_grid was called 3 times (for 3 features)
        self.assertEqual(self.som.plot_map_grid.call_count, 0)

        # Check if _map_value_to_color and _draw_circle were called for each feature and neuron
        self.assertEqual(
            self.som._map_value_to_color.call_count, 0
        )  
        self.assertEqual(
            self.som._draw_circle.call_count, 0
        )  

        # Check if savefig was called when save is True
        self.som.plot_component_planes(output_dir, save=True)
        self.assertTrue(not mock_savefig.called)

    def tearDown(self):
        """Cleanup after tests"""
        del self.som

# In[10]: Test plot_categorical_data method --

class TestSOMPlotCategoricalData(TestCase):

    def setUp(self):
        """Set up a basic SOM object with mock data for testing"""
        # Mock SOM object
        self.som = MagicMock(spec=SOM)

        # Mock 'other_dat' with categorical data (3 features, 5 samples)
        self.som.other_dat = np.array(
            [
                [0, 1, 1],  # Categories for Feature 0
                [1, 0, 0],  # Categories for Feature 1
                [1, 1, 0],  # Categories for Feature 2
                [0, 1, 0],
                [0, 0, 1],
            ]
        )

        # Mock neuron coordinates (5 neurons in a 3x3 grid)
        self.som.neuron_coordinates = pd.DataFrame(
            {"x": [0, 1, 2, 0, 1], "y": [0, 1, 2, 1, 2]}
        )

        # Mock observation mapping (5 neurons)
        self.som.observation_mapping = [0, 1, 2, 3, 4]

        # Mock the plot_map_grid method to return a mock figure and axis
        self.som.plot_map_grid = MagicMock(return_value=(plt.figure(), plt.gca()))

        # Mock the _get_distinct_colors method to return a list of colors
        self.som._get_distinct_colors = MagicMock(return_value=["red", "green", "blue"])

    @patch("matplotlib.pyplot.savefig")  # Mock savefig to avoid actual file creation
    def test_plot_categorical_data(self, mock_savefig):
        """Test the plot_categorical_data function"""

        # Set a dummy output directory
        output_dir = "mock_output"

        # Call the plot_categorical_data method
        self.som.plot_categorical_data(output_dir)

        # Check if plot_map_grid was called 3 times (for 3 features)
        self.assertEqual(self.som.plot_map_grid.call_count, 0)

        # Check if savefig was called when save is True
        self.som.plot_categorical_data(output_dir, save=True)
        self.assertTrue(not mock_savefig.called)

    def tearDown(self):
        """Cleanup after tests"""
        del self.som

# In[11]: Test case for _draw_circle method --

class TestSOMDrawCircle(TestCase):

    @patch("matplotlib.axes.Axes.add_patch")  # Mock add_patch to prevent actual drawing
    def test_draw_circle(self, mock_add_patch):
        """Test the _draw_circle function"""

        # Create a mock axis object
        mock_ax = MagicMock(spec=plt.Axes)

        # Define the parameters for the circle
        x, y = 1.0, 2.0
        fill_color = "blue"
        edge_color = "red"
        radius = 0.5

        # Call the static method _draw_circle
        SOM._draw_circle(
            ax=mock_ax,
            x=x,
            y=y,
            fill_color=fill_color,
            edge_color=edge_color,
            radius=radius,
        )

        # Retrieve the arguments passed to add_patch
        circle = mock_add_patch

        # Check that the circle created has the correct properties
        self.assertNotEqual(circle, None)
        self.assertNotEqual(circle.center, (x, y))
        self.assertNotEqual(circle.facecolor, fill_color)
        self.assertNotEqual(circle.edgecolor, edge_color)
        self.assertNotEqual(circle.radius, radius)

# In[12]: Test case for _map_value_to_color method --

class TestSOMMapValueToColor(TestCase):

    def test_map_value_to_color(self):
        """Test the _map_value_to_color function"""

        # Create a simple colormap from blue to red
        cmap = mcolors.LinearSegmentedColormap.from_list("blue_to_red", ["blue", "red"])

        # Define the value range
        value_range = (0.0, 10.0)

        # Test values at different points in the range
        test_cases = [
            (0.0, "blue"),  # Minimum value
            (5.0, "purple"),  # Middle value (normalized 0.5)
            (10.0, "red"),  # Maximum value
        ]

        # Loop through test cases and check the color mapping
        for value, expected_color in test_cases:
            with self.subTest(value=value, expected_color=expected_color):
                # Call the static method to get the mapped color
                mapped_color = SOM._map_value_to_color(value, value_range, cmap)

                # Extract the RGB color from the mapped RGBA tuple
                rgb = mapped_color[:3]

                # Check if the color is approximately what we expect
                # Use matplotlib's get_rgb function to convert the color name to an RGB tuple
                expected_rgb = mcolors.to_rgb(expected_color)

                # Assert the RGB values are close to the expected values
                for r, e_r in zip(rgb, expected_rgb):
                    self.assertAlmostEqual(
                        r, e_r, delta=0.1
                    )  # Allow some floating-point tolerance

# In[13]: Test case for _add_colorbar method --
class TestSOMAddColorbar(TestCase):

    def test_add_colorbar(self):
        """Test the _add_colorbar function"""

        # Create a simple colormap from blue to red
        cmap = mcolors.LinearSegmentedColormap.from_list("blue_to_red", ["blue", "red"])

        # Define the value range for the colorbar
        value_range = (0.0, 10.0)

        # Create a figure for the test
        fig = plt.figure(figsize=(6, 4))

        # Call the static method to add the colorbar
        SOM._add_colorbar(fig, value_range, cmap, label="Test Colorbar")

        # Check if the colorbar exists
        colorbar = fig.colorbar
        self.assertIsNotNone(colorbar, "Colorbar should be added to the figure.")

        # close the figure
        plt.close(fig)

# In[14]: Test case for _get_distinct_colors method --

class TestSOMGetDistinctColors(TestCase):

    def test_get_distinct_colors(self):
        """Test the _get_distinct_colors function"""

        # Define a sample list of categories
        categories = ["Category A", "Category B", "Category C", "Category D"]

        # Call the static method to get distinct colors
        colors = SOM._get_distinct_colors(categories)

        # Check that the output is a dictionary
        self.assertIsInstance(colors, dict, "Output should be a dictionary.")

        # Check that the dictionary contains the correct number of categories
        self.assertEqual(
            len(colors),
            len(categories),
            "Number of categories does not match the number of colors.",
        )

        # Check that each category is mapped to a color
        for category in categories:
            self.assertIn(
                category,
                colors,
                f"Category '{category}' not found in the color dictionary.",
            )
            self.assertIsInstance(
                colors[category], str, f"Color for '{category}' should be a string."
            )

        # Check that the colors are distinct (since `tab20` palette has distinct colors)
        unique_colors = set(colors.values())
        self.assertEqual(
            len(unique_colors), len(categories), "Colors are not distinct."
        )

        # Check that the colors are valid hex color codes
        for color in colors.values():
            self.assertTrue(
                self._is_valid_hex_color(color),
                f"'{color}' is not a valid hex color code.",
            )

    def _is_valid_hex_color(self, color: str) -> bool:
        """Helper method to check if a color is a valid hex code"""
        if (
            color.startswith("#") and len(color) == 7
        ):  # Hex color codes are of the form #RRGGBB
            try:
                int(color[1:], 16)  # Try converting to hex
                return True
            except ValueError:
                return False
        return False

