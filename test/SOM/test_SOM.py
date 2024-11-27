# Standard imports
import os
import sys

# Third party imports
import numpy as np
import pandas as pd
import pytest

# Local imports
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from SOM.utils.som_utils import SOM


class TestSOM:
    """
    Test suite for the SOM class.
    """

    @pytest.fixture
    def sample_data(self):
        """
        Fixture to provide sample training and categorical data with 100 observations.
        """
        # Generate synthetic training data with 100 observations and 3 features
        np.random.seed(42)  # For reproducibility
        train_data = pd.DataFrame({
            "feature1": np.random.randint(1, 10, size=100),
            "feature2": np.random.randint(2, 10, size=100),
            "feature3": np.random.randint(3, 10, size=100),
        })

        # Generate corresponding categorical data
        other_data = pd.DataFrame({
            "category": np.random.choice(["A", "B"], size=100)
        })

        return train_data, other_data


    @pytest.fixture
    def som_instance(self, sample_data):
        """
        Fixture to provide a basic SOM instance.
        """
        train_data, other_data = sample_data
        som = SOM(
            train_dat=train_data,
            other_dat=other_data,
            scale_method="zscore",
            x_dim=3,
            y_dim=2,
            topology="rectangular",
            neighborhood_fnc="gaussian",
            epochs=10
        )
        return som


    def test_som_initialization(self, som_instance):
        """
        Test that the SOM instance initializes properly.
        """
        assert som_instance.xdim == 3
        assert som_instance.ydim == 2
        assert som_instance.topology == "hexagonal"
        assert som_instance.neighborhood_fnc == "gaussian"
        assert som_instance.epochs == 10
        assert som_instance.map is None
        assert som_instance.weights is None
        assert som_instance.weights_scaled is None
        assert som_instance.observation_mapping is None
        assert som_instance.neuron_coordinates is None


    def test_som_train_map(self, som_instance):
        """
        Test SOM training.
        """
        som_instance.train_map()
        assert som_instance.map is not None
        assert som_instance.weights is not None
        assert som_instance.weights_scaled is not None
        assert som_instance.observation_mapping is not None
        assert som_instance.neuron_coordinates is not None


    def test_non_numeric_training_data(self):
        """
        Test that TypeError is returned for non-numeric training data.
        """
        train_dat = pd.DataFrame({
            "feature1": [1, "A", 9, 4, 7, 1],
            "feature2": [1, 8, 2, 2, 3, 7],
            "feature3": [4, 6, 7, 9, "B", 4]
        })
        other_dat = pd.DataFrame({
            "category": ["A", "B", "A", "A", "A", "B"]
        })

        with pytest.raises(
            expected_exception=TypeError,
            match="train_dat must contain only numeric values."
        ):
            SOM(
                train_dat=train_dat,
                other_dat=other_dat,
                scale_method="zscore",
                x_dim=3,
                y_dim=2,
                topology="hexagonal",
                neighborhood_fnc="gaussian",
                epochs=10
            )


    def test_empty_training_data(self):
        """
        Test that ValueError is returned when training data is empty.
        """
        train_dat = pd.DataFrame()
        other_dat = pd.DataFrame()

        with pytest.raises(
            expected_exception=ValueError,
            match="Training data must have at least two features."
        ):
            SOM(
                train_dat=train_dat,
                other_dat=other_dat,
                scale_method="zscore",
                x_dim=3,
                y_dim=2,
                topology="hexagonal",
                neighborhood_fnc="gaussian",
                epochs=10
            )


    def test_unequal_data_sizes(self, sample_data):
        """
        Test that a ValueError is raised when the number of rows in `train_dat` and `other_dat` do
        not match.
        """
        train_data, other_data = sample_data

        # Modify `other_data` to have a different number of rows
        other_data_mod = other_data.iloc[:-1]  # Drop the last row

        with pytest.raises(
            expected_exception=ValueError,
            match="The number of rows in `train_dat` must match the number of rows in `other_dat`."
        ):
            SOM(
                train_dat=train_data,
                other_dat=other_data_mod,
                scale_method="zscore",
                x_dim=5,
                y_dim=5,
                topology="hexagonal",
                neighborhood_fnc="gaussian",
                epochs=100
            )


    def test_invalid_scaling_method(self, sample_data):
        """
        Test that ValueError is returned for invalid scaling method.
        """
        train_data, other_data = sample_data
        with pytest.raises(
            expected_exception=ValueError,
            match="scale_method must be 'zscore' or 'minmax'."
        ):
            SOM(
                train_dat=train_data,
                other_dat=other_data,
                scale_method="invalid",
                x_dim=3,
                y_dim=2,
                topology="hexagonal",
                neighborhood_fnc="gaussian",
                epochs=10
            )


    def test_invalid_x_dimension(self, sample_data):
        """
        Testthat ValueError is returned invalid x-dimension (negative value).
        """
        train_data, other_data = sample_data
        with pytest.raises(
            expected_exception=ValueError,
            match="x_dim must be a positive integer."
        ):
            SOM(
                train_dat=train_data,
                other_dat=other_data,
                scale_method="zscore",
                x_dim=-1,
                y_dim=2,
                topology="hexagonal",
                neighborhood_fnc="gaussian",
                epochs=10
            )


    def test_invalid_y_dimension(self, sample_data):
        """
        Test that ValueError is returned for invalid y-dimension (non-integer).
        """
        train_data, other_data = sample_data
        with pytest.raises(
            expected_exception=ValueError,
            match="y_dim must be a positive integer."
        ):
            SOM(
                train_dat=train_data,
                other_dat=other_data,
                scale_method="zscore",
                x_dim=2,
                y_dim=2.4,
                topology="hexagonal",
                neighborhood_fnc="gaussian",
                epochs=10
            )


    def test_calculate_percent_variance_explained(self, som_instance):
        """
        Test Percent Variance Explained calculation.
        """
        som_instance.train_map()
        pve = som_instance.calculate_percent_variance_explained()

        # PVE should be a value between 0 and 100
        assert 0 <= pve <= 100


    def test_calculate_topographic_error(self, som_instance):
        """
        Test Topographic Error calculation.
        """
        som_instance.train_map()
        topo_error = som_instance.calculate_topographic_error()

        # Topographic error should be a value between 0 and 1
        assert 0 <= topo_error <= 1


    def test_zscore_scaling(self, som_instance):
        """
        Test that scaling and unscaling with Z-score normalization returns the original data.
        """
        # Use the original training data from som_instance
        original_data = som_instance.train_dat

        # Perform Z-score scaling using SOM instance
        scaled_data, scaling_factors = som_instance._zscore_scale()

        # Set scaling factors explicitly to ensure unscale uses the same
        som_instance._scaling_factors = scaling_factors

        # Reverse the scaling
        unscaled_data = som_instance._zscore_unscale(scaled_data)

        # Assert that the unscaled data matches the original data
        np.testing.assert_almost_equal(
            unscaled_data, original_data, decimal=5
        )


    def test_minmax_scaling(self, som_instance):
        """
        Test that scaling and unscaling with min-max normalization returns the original data.
        """
        # Use the original training data from som_instance
        original_data = som_instance.train_dat

        # Perform min-max scaling using SOM instance
        scaled_data, scaling_factors = som_instance._minmax_scale()

        # Set scaling factors explicitly to ensure unscale uses the same
        som_instance._scaling_factors = scaling_factors

        # Reverse the scaling
        unscaled_data = som_instance._minmax_unscale(scaled_data)

        # Assert that the unscaled data matches the original data
        np.testing.assert_almost_equal(
            unscaled_data, original_data, decimal=5
        )


    def test_plot_component_planes(self, som_instance, tmp_path):
        """
        Test plotting component planes.
        """
        som_instance.train_map()
        output_dir = tmp_path / "plots"
        som_instance.plot_component_planes(output_dir=str(output_dir))

        # The number of component planes generated should equal the number of features in the
        # training data
        assert len(list(output_dir.iterdir())) == 3


    def test_plot_categorical_data(self, som_instance, tmp_path):
        """
        Test plotting categorical data.
        """
        som_instance.train_map()
        output_dir = tmp_path / "plots"
        som_instance.plot_categorical_data(output_dir=str(output_dir))

        # The number of categorical plots generated should equal the number of features in the
        # other/categorical data
        assert len(list(output_dir.iterdir())) == 1
