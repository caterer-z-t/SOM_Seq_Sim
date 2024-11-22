import pytest
import pandas as pd
from unittest.mock import patch
import os
import numpy as np
import yaml
import unittest

from Seq_Sim.utils.seq_sim_utils import (
    save_data_to_csv, 
    create_directory_if_not_exists,
    validate_arguments,
    case_when,
    load_config,
    calculate_diff_expression,
    generate_diff_cell_types,
    identify_diff_clusters,
    create_new_rows_for_cell_type,
    generate_rare_cell_counts,
    generate_major_cell_counts,
    apply_noise,
    encode_categorical_columns,
    set_random_seed,
    combine_features,
    generate_variance,
    generate_cluster_features,
    add_rows_for_diff_cells
)

# Test class using unittest
class TestAddRowsForDiffCells(unittest.TestCase):

    """Test the add_rows_for_diff_cells function.
    """

    def setUp(self):
        """Prepare mock data for the tests.
        """
        # Prepare mock data for the tests
        self.celltype_df = pd.DataFrame({
            "subject_id": [1, 2, 3],
            "cell_type": ["A", "B", "A"]
        })

        self.dummy_data = pd.DataFrame({
            "subject_id": [1, 2, 3],
            "disease": [1, 0, 1]  # Subjects 1 and 3 have the disease
        })

    def test_add_rows_for_diff_cells_no_disease(self):
        """Test the function with no subjects having the disease.
        """
        # Test with diff = 2, and no subjects with disease
        dummy_data_no_disease = pd.DataFrame({
            "subject_id": [1, 2, 3],
            "disease": [0, 0, 0]  # No subjects with disease
        })
        diff = 2
        cell_type = "C"
        
        result = add_rows_for_diff_cells(self.celltype_df, dummy_data_no_disease, diff, cell_type)

        # Expected: No new rows should be added
        pd.testing.assert_frame_equal(result, self.celltype_df)

# Test class for generate_cluster_features
class TestGenerateClusterFeatures(unittest.TestCase):
    """Test suite for the generate_cluster_features function.

    Args:
        unittest (TestCase): The base class for test cases.
    """

    def test_generate_cluster_features_normal(self):
        """Test generate_cluster_features with normal data
    """
        cell_clusters = np.array([0, 1, 2, 1, 0, 2])  # Sample encoded clusters
        unique_clusters = 3  # Number of unique clusters

        # Expected variance and generated cluster features
        features, variance = generate_cluster_features(cell_clusters, unique_clusters)

        # Check that the features and variance are arrays and have the same length as cell_clusters
        self.assertEqual(len(features), len(cell_clusters))
        
        # Check that the variance is a positive value
        self.assertGreater(variance, 0)

    def test_generate_cluster_features_empty(self):
        """Test generate_cluster_features with an empty array"""
        cell_clusters = np.array([])  # No clusters
        unique_clusters = 0  # No unique clusters

        features, _ = generate_cluster_features(cell_clusters, unique_clusters)

        # Check that it handles empty arrays gracefully
        self.assertEqual(len(features), 0)
        

    def test_generate_cluster_features_single_cluster(self):
        """Test generate_cluster_features with a single cluster"""
        cell_clusters = np.array([0, 0, 0, 0, 0])  # All cells in the same cluster
        unique_clusters = 1  # Only one unique cluster

        features, variance = generate_cluster_features(cell_clusters, unique_clusters)

        # Check that variance is calculated and features are generated
        self.assertEqual(len(features), len(cell_clusters))
        self.assertGreater(variance, 0)

# Test class for generate_variance
class TestGenerateVariance(unittest.TestCase):
    
    def test_generate_variance_normal_data(self):
        """Test generate_variance with normal data"""
        data = np.array([1, 2, 3, 4, 5])
        result = generate_variance(data)
        expected = np.mean(1 / (data + 1e-8))
        self.assertAlmostEqual(result, expected, places=6)

    def test_generate_variance_large_values(self):
        """Test generate_variance with very large values"""
        data = np.array([1e6, 1e7, 1e8])
        result = generate_variance(data)
        expected = np.mean(1 / (data + 1e-8))
        self.assertAlmostEqual(result, expected, places=6)

    def test_generate_variance_single_value(self):
        """Test generate_variance with a single value"""
        data = np.array([1])
        result = generate_variance(data)
        expected = 1 / (data + 1e-8)  # Variance is just the inverse of the value
        self.assertAlmostEqual(result, expected, places=6)

# Test class for combine_features
class TestCombineFeatures(unittest.TestCase):

    def test_empty_features(self):
        """Test the case when features list is empty"""
        features = []
        n_features = 0
        
        # Combine the empty features
        features_df = combine_features(features, n_features)
        
        # Assert the DataFrame is empty
        self.assertTrue(features_df.empty, "DataFrame should be empty for empty features.")

    def test_single_feature(self):
        """Test combining a single feature into a DataFrame"""
        features = [[1]]
        n_features = 1
        
        # Combine the features
        features_df = combine_features(features, n_features)
        
        # Assert the shape is correct
        self.assertEqual(features_df.shape, (1, 1), "Shape of the DataFrame is incorrect for a single feature.")
        
        # Assert the column name is correct
        self.assertListEqual(list(features_df.columns), ["feature_1"], "Column name is incorrect for a single feature.")

# Test class for set_random_seed
class TestSetRandomSeed(unittest.TestCase):

    def test_set_random_seed(self):
        """Test setting random seed for reproducibility"""
        # Set the random seed using specific values for x and seed
        x = 5
        seed = 10
        
        # Save the result of a random number generated before setting the seed
        np.random.seed(None)  # Reset the seed before the test
        random_before = np.random.random()
        
        # Set the random seed using the function
        set_random_seed(x, seed)
        
        # Save the result of a random number generated after setting the seed
        random_after = np.random.random()
        
        # Check that the generated random values are consistent after setting the seed
        set_random_seed(x, seed)  # Set the same seed again
        random_after_repeat = np.random.random()
        
        # Assert that the same seed produces the same random number
        self.assertEqual(random_after, random_after_repeat, "Random values are not reproducible.")
        
        # Assert that the random number changes when the seed is not set
        self.assertNotEqual(random_before, random_after, "Random values should change before and after setting the seed.")

    def test_edge_case_seed(self):
        """Test edge case when the seed is zero"""
        x = 5
        seed = 0
        
        # Set random seed with zero seed
        set_random_seed(x, seed)
        
        # Generate a random value and check its consistency
        np.random.seed(x * seed)  # Reset seed for comparison
        random_value = np.random.random()
        
        # Check that the random value is reproducible
        set_random_seed(x, seed)
        random_value_repeat = np.random.random()
        
        self.assertEqual(random_value, random_value_repeat, "Random values are not reproducible with zero seed.")
        
    def test_large_seed(self):
        """Test edge case with large seed value"""
        x = 5
        seed = 1_000_000
        
        # Set random seed with a large seed value
        set_random_seed(x, seed)
        
        # Generate a random value and check its consistency
        np.random.seed(x * seed)  # Reset seed for comparison
        random_value = np.random.random()
        
        # Check that the random value is reproducible
        set_random_seed(x, seed)
        random_value_repeat = np.random.random()
        
        self.assertEqual(random_value, random_value_repeat, "Random values are not reproducible with large seed.")

# Test class for encode_categorical_columns
class TestEncodeCategoricalColumns(unittest.TestCase):

    def setUp(self):
        """Setup dummy data for testing"""
        self.data = pd.DataFrame({
            'cluster': ['A', 'B', 'A', 'C', 'B'],
            'disease': ['Disease1', 'Disease2', 'Disease1', 'Disease3', 'Disease2'],
            'individual': ['Ind1', 'Ind2', 'Ind1', 'Ind3', 'Ind2']
        })

    def test_encode_categorical_columns(self):
        """Test encoding of categorical columns"""
        cell_diseases, cell_individual, cell_clusters = encode_categorical_columns(
            self.data, 'cluster', 'disease', 'individual'
        )
        
        # Check that the result is a tuple of three arrays
        self.assertEqual(len(cell_diseases), len(self.data))
        self.assertEqual(len(cell_individual), len(self.data))
        self.assertEqual(len(cell_clusters), len(self.data))
        
        # Check that all values are numeric (encoded)
        self.assertTrue(np.issubdtype(cell_diseases.dtype, np.integer))
        self.assertTrue(np.issubdtype(cell_individual.dtype, np.integer))
        self.assertTrue(np.issubdtype(cell_clusters.dtype, np.integer))
        
        # Check correct encoding for the categorical columns
        expected_diseases = [0, 1, 0, 2, 1]  # Disease1 -> 0, Disease2 -> 1, Disease3 -> 2
        expected_individuals = [0, 1, 0, 2, 1]  # Ind1 -> 0, Ind2 -> 1, Ind3 -> 2
        expected_clusters = [0, 1, 0, 2, 1]  # A -> 0, B -> 1, C -> 2
        
        np.testing.assert_array_equal(cell_diseases, expected_diseases)
        np.testing.assert_array_equal(cell_individual, expected_individuals)
        np.testing.assert_array_equal(cell_clusters, expected_clusters)

    def test_empty_data(self):
        """Test empty DataFrame"""
        empty_data = pd.DataFrame(columns=['cluster', 'disease', 'individual'])
        cell_diseases, cell_individual, cell_clusters = encode_categorical_columns(
            empty_data, 'cluster', 'disease', 'individual'
        )
        
        # Ensure that empty data results in empty arrays
        self.assertEqual(len(cell_diseases), 0)
        self.assertEqual(len(cell_individual), 0)
        self.assertEqual(len(cell_clusters), 0)

    def test_missing_column(self):
        """Test missing column"""
        with self.assertRaises(KeyError):
            encode_categorical_columns(self.data, 'cluster', 'missing_column', 'individual')

# Test class for apply_noise
class TestApplyNoise(unittest.TestCase):

    def test_basic_case(self):
        """Test with typical values for the input parameters."""
        np.random.seed(42)  # Set seed for reproducibility
        result = apply_noise(5, [0.1, 0.2, 0.3], 10, 0.5, 0.1)
        
        # Check that the result is a numpy array of the correct shape
        self.assertEqual(result.shape, (10,))
        
        # Check that the final feature contains noise (non-zero variance)
        self.assertFalse(np.all(result == 5 * 0.5))  # This would be the result without noise

    def test_zero_cells(self):
        """Test with zero cells, expecting an empty result."""
        result = apply_noise(5, [0.1, 0.2, 0.3], 0, 0.5, 0.1)
        self.assertEqual(result.shape, (0,))

    def test_edge_case_negative_cluster_ratio(self):
        """Test with negative cluster ratio, expecting a valid result."""
        np.random.seed(42)
        result = apply_noise(5, [0.1, 0.2, 0.3], 10, -0.5, 0.1)
        
        # Check that the result has the correct shape
        self.assertEqual(result.shape, (10,))

    def test_edge_case_large_cluster_ratio(self):
        """Test with a very large cluster ratio."""
        np.random.seed(42)
        result = apply_noise(5, [0.1, 0.2, 0.3], 10, 1000, 0.1)
        
        # Check that the result has the correct shape
        self.assertEqual(result.shape, (10,))
        
        # Check that the result contains a large signal component due to the large cluster ratio
        self.assertTrue(np.all(result > 500))

    def test_large_number_of_cells(self):
        """Test with a large number of cells."""
        np.random.seed(42)
        result = apply_noise(5, [0.1, 0.2, 0.3], 10000, 0.5, 0.1)
        
        # Check that the result is a numpy array with 10000 values
        self.assertEqual(result.shape, (10000,))

# Test class for generate_major_cell_counts
class TestGenerateMajorCellCounts(unittest.TestCase):

    def test_basic_case(self):
        """Test with typical values for the input parameters."""
        np.random.seed(42)  # Set seed for reproducibility
        result = generate_major_cell_counts(1000, 0.1, 5)
        # Check that the result is a list of 5 values
        self.assertEqual(len(result), 5)
        # Check that values are within the expected range
        for val in result:
            self.assertGreaterEqual(val, 1000 * (1 - 0.1))
            self.assertLessEqual(val, 1000 * (1 + 0.1))

    def test_large_number_of_major_cell_types(self):
        """Test with a large number of major cell types."""
        np.random.seed(42)
        result = generate_major_cell_counts(1000, 0.1, 1000)
        self.assertEqual(len(result), 1000)
        for val in result:
            self.assertGreaterEqual(val, 1000 * (1 - 0.1))
            self.assertLessEqual(val, 1000 * (1 + 0.1))

    def test_edge_case_standard_deviation_zero(self):
        """Test when the standard deviation of cell types is zero."""
        np.random.seed(42)
        result = generate_major_cell_counts(1000, 0.0, 5)
        # Check if all the values are equal because standard deviation is zero
        self.assertTrue(np.all(result == 1000))

    def test_non_integer_values(self):
        """Test with non-integer values for n_cells, expecting rounding to the nearest integer."""
        np.random.seed(42)
        result = generate_major_cell_counts(1000.5, 0.1, 5)
        self.assertEqual(len(result), 5)
        # Check if values are rounded and within the expected range
        for val in result:
            self.assertGreaterEqual(val, 1000 * (1 - 0.1))
            self.assertLessEqual(val, 1000 * (1 + 0.1))

# Test class for generate_rare_cell_counts
class TestGenerateRareCellCounts(unittest.TestCase):

    def test_basic_case(self):
        """Test with typical values for the input parameters."""
        np.random.seed(42)  # Set seed for reproducibility
        result = generate_rare_cell_counts(1000, 0.1, 0.05, 5)
        # Check that the result is a list of 5 values
        self.assertEqual(len(result), 5)
        # Check that values are within the expected range
        for val in result:
            self.assertGreaterEqual(val, 1000 * 0.05 * (1 - 0.1))
            self.assertLessEqual(val, 1000 * 0.05 * (1 + 0.1))

    def test_large_number_of_minor_cell_types(self):
        """Test with a large number of minor cell types."""
        np.random.seed(42)
        result = generate_rare_cell_counts(1000, 0.1, 0.05, 1000)
        self.assertEqual(len(result), 1000)
        for val in result:
            self.assertGreaterEqual(val, 1000 * 0.05 * (1 - 0.1))
            self.assertLessEqual(val, 1000 * 0.05 * (1 + 0.1))

    def test_edge_case_standard_deviation_zero(self):
        """Test when the standard deviation of cell types is zero."""
        np.random.seed(42)
        result = generate_rare_cell_counts(1000, 0.0, 0.05, 5)
        # Check if all the values are equal because standard deviation is zero
        self.assertTrue(np.all(result == 1000 * 0.05))

    def test_edge_case_relative_abundance_zero(self):
        """Test when relative abundance is zero."""
        np.random.seed(42)
        result = generate_rare_cell_counts(1000, 0.1, 0.0, 5)
        # All cell counts should be zero since relative abundance is zero
        self.assertTrue(np.all(result == 0))

    def test_non_integer_values(self):
        """Test with non-integer values for n_cells, expecting rounding to the nearest integer."""
        np.random.seed(42)
        result = generate_rare_cell_counts(1000.5, 0.1, 0.05, 5)
        self.assertEqual(len(result), 5)
        # Check if values are rounded and within the expected range
        for val in result:
            self.assertGreaterEqual(val, 1000 * 0.05 * (1 - 0.1))
            self.assertLessEqual(val, 1000 * 0.05 * (1 + 0.1))

# Test class for create_new_rows_for_cell_type
class TestCreateNewRowsForCellType(unittest.TestCase):

    def test_basic_case(self):
        """Test with a typical case for subject ID, cell type, and count."""
        result = create_new_rows_for_cell_type("S1", "T-cell", 3)
        expected = pd.DataFrame({
            "cell_type": ["T-cell", "T-cell", "T-cell"],
            "subject_id": ["S1", "S1", "S1"]
        })
        pd.testing.assert_frame_equal(result, expected)

    def test_large_count(self):
        """Test with a large count to check performance."""
        result = create_new_rows_for_cell_type("S1", "T-cell", 1000)
        expected = pd.DataFrame({
            "cell_type": ["T-cell"] * 1000,
            "subject_id": ["S1"] * 1000
        })
        pd.testing.assert_frame_equal(result, expected)

    def test_non_integer_count(self):
        """Test with a non-integer count, expecting the count to be truncated."""
        result = create_new_rows_for_cell_type("S1", "T-cell", 3.7)
        expected = pd.DataFrame({
            "cell_type": ["T-cell", "T-cell", "T-cell"],
            "subject_id": ["S1", "S1", "S1"]
        })
        pd.testing.assert_frame_equal(result, expected)

    def test_different_cell_type(self):
        """Test with a different cell type."""
        result = create_new_rows_for_cell_type("S1", "B-cell", 2)
        expected = pd.DataFrame({
            "cell_type": ["B-cell", "B-cell"],
            "subject_id": ["S1", "S1"]
        })
        pd.testing.assert_frame_equal(result, expected)

# Test class for identify_diff_clusters
class TestIdentifyDiffClusters(unittest.TestCase):

    def test_basic_case(self):
        """Test with typical values for major and minor differentially expressed cell types."""
        result = identify_diff_clusters(3, 2, 10)
        expected = [0, 1, 2, 8, 9]  # First 3 for major, last 2 for minor
        self.assertEqual(result, expected)

    def test_no_major_diff_celltypes(self):
        """Test with no major differentially expressed cell types."""
        result = identify_diff_clusters(0, 2, 10)
        expected = [8, 9]  # Only minor cell types
        self.assertEqual(result, expected)

    def test_no_minor_diff_celltypes(self):
        """Test with no minor differentially expressed cell types."""
        result = identify_diff_clusters(3, 0, 10)
        expected = [0, 1, 2]  # Only major cell types
        self.assertEqual(result, expected)

    def test_no_differential_cell_types(self):
        """Test when no differential cell types are specified."""
        result = identify_diff_clusters(0, 0, 10)
        expected = []  # No differentially expressed clusters
        self.assertEqual(result, expected)

    def test_large_input(self):
        """Test with larger values for major and minor differentially expressed cell types."""
        result = identify_diff_clusters(5, 5, 20)
        expected = list(range(5)) + list(range(15, 20))  # First 5 for major, last 5 for minor
        self.assertEqual(result, expected)

# Test class for generate_diff_cell_types
class TestGenerateDiffCellTypes(unittest.TestCase):

    def test_basic_case(self):
        """Test with a standard list of cluster indices."""
        diff_clusters = [0, 1, 2]
        result = generate_diff_cell_types(diff_clusters)
        expected = ['A', 'B', 'C']
        self.assertEqual(result, expected)

    def test_empty_list(self):
        """Test with an empty list of cluster indices."""
        diff_clusters = []
        result = generate_diff_cell_types(diff_clusters)
        expected = []
        self.assertEqual(result, expected)

    def test_large_input(self):
        """Test with a large list of indices."""
        diff_clusters = list(range(50))
        result = generate_diff_cell_types(diff_clusters)
        # Expecting 'A' to 'AX' as the first 50 letters
        expected = [chr(65 + i) for i in diff_clusters]
        self.assertEqual(result, expected)

# Test class for calculate_diff_expression
class TestCalculateDiffExpression(unittest.TestCase):
    def test_zero_fold_change(self):
        """Test with a fold change of zero."""
        abundance = np.array([0.1, 0.2, 0.15])
        n_cells = 100
        fc_interact = 0.0
        result = calculate_diff_expression(abundance, n_cells, fc_interact)
        expected = 0  # No additional rows when fold change is zero
        self.assertEqual(result, expected)

    def test_negative_fold_change(self):
        """Test with a negative fold change."""
        abundance = np.array([0.3, 0.4, 0.35])
        n_cells = 100
        fc_interact = -0.2
        result = calculate_diff_expression(abundance, n_cells, fc_interact)
        expected = -7  # Expected from manual calculation
        self.assertEqual(result, expected)

    def test_single_value_abundance(self):
        """Test with a single value in abundance."""
        abundance = np.array([0.25])
        n_cells = 100
        fc_interact = 0.3
        result = calculate_diff_expression(abundance, n_cells, fc_interact)
        expected = 8  # Expected from manual calculation
        self.assertEqual(result, expected)

# Test class for the load_config function
class TestLoadConfig(unittest.TestCase):
    def setUp(self):
        """Create temporary YAML files for testing."""
        self.valid_config_path = "valid_config.yaml"
        self.empty_config_path = "empty_config.yaml"
        self.invalid_config_path = "invalid_config.yaml"

        # Create a valid YAML configuration file
        valid_config = {
            "setting1": "value1",
            "setting2": 42,
            "setting3": [1, 2, 3]
        }
        with open(self.valid_config_path, "w") as file:
            yaml.dump(valid_config, file)

        # Create an empty YAML file
        with open(self.empty_config_path, "w") as file:
            file.write("")

        # Create an invalid YAML file
        with open(self.invalid_config_path, "w") as file:
            file.write("key: value: another_value")

    def tearDown(self):
        """Remove temporary files after testing."""
        os.remove(self.valid_config_path)
        os.remove(self.empty_config_path)
        os.remove(self.invalid_config_path)

    def test_valid_config(self):
        """Test loading a valid YAML configuration."""
        expected_config = {
            "setting1": "value1",
            "setting2": 42,
            "setting3": [1, 2, 3]
        }
        result = load_config(self.valid_config_path)
        self.assertEqual(result, expected_config)

    def test_empty_config(self):
        """Test loading an empty YAML configuration."""
        result = load_config(self.empty_config_path)
        self.assertIsNone(result)

    def test_invalid_config(self):
        """Test loading an invalid YAML configuration."""
        with self.assertRaises(yaml.YAMLError):
            load_config(self.invalid_config_path)

    def test_file_not_found(self):
        """Test loading a non-existent YAML file."""
        with self.assertRaises(FileNotFoundError):
            load_config("non_existent.yaml")

# Test class for the case_when function
class TestCaseWhen(unittest.TestCase):
    def test_single_condition_true(self):
        """Test with a single true condition."""
        conditions = [(True, "Matched")]
        default = "Default"
        result = case_when(conditions, default)
        self.assertEqual(result, "Matched")

    def test_multiple_conditions_first_true(self):
        """Test with multiple conditions where the first is true."""
        conditions = [(True, "First"), (True, "Second"), (False, "Third")]
        default = "Default"
        result = case_when(conditions, default)
        self.assertEqual(result, "First")

    def test_multiple_conditions_second_true(self):
        """Test with multiple conditions where the second is true."""
        conditions = [(False, "First"), (True, "Second"), (False, "Third")]
        default = "Default"
        result = case_when(conditions, default)
        self.assertEqual(result, "Second")

    def test_all_conditions_false(self):
        """Test with all conditions false."""
        conditions = [(False, "First"), (False, "Second"), (False, "Third")]
        default = "Default"
        result = case_when(conditions, default)
        self.assertEqual(result, "Default")

    def test_empty_conditions_list(self):
        """Test with an empty conditions list."""
        conditions = []
        default = "Default"
        result = case_when(conditions, default)
        self.assertEqual(result, "Default")

    def test_mixed_conditions(self):
        """Test with mixed conditions including complex evaluations."""
        conditions = [(2 > 1, "Greater"), (3 < 1, "Lesser"), (4 == 4, "Equal")]
        default = "Default"
        result = case_when(conditions, default)
        self.assertEqual(result, "Greater")

# Test class for the validate_arguments function
class TestValidateArguments(unittest.TestCase):
    def test_valid_arguments(self):
        """Test with valid arguments."""
        args = ["my_script.py", "10", "1.5", "config.json"]
        expected = (10, 1.5, "config.json")
        result = validate_arguments(args)
        self.assertEqual(result, expected)

    def test_missing_arguments(self):
        """Test with missing arguments."""
        args = ["my_script.py", "10"]
        with self.assertRaises(ValueError) as context:
            validate_arguments(args)
        self.assertIn("Usage: python my_script.py", str(context.exception))

    def test_invalid_num_samples(self):
        """Test with invalid num_samples."""
        args = ["my_script.py", "ten", "1.5", "config.json"]
        with self.assertRaises(ValueError) as context:
            validate_arguments(args)
        self.assertIn("Invalid argument type", str(context.exception))

    def test_invalid_fold_change(self):
        """Test with invalid fold_change."""
        args = ["my_script.py", "10", "one.point.five", "config.json"]
        with self.assertRaises(ValueError) as context:
            validate_arguments(args)
        self.assertIn("Invalid argument type", str(context.exception))

    def test_extra_arguments(self):
        """Test with extra arguments."""
        args = ["my_script.py", "10", "1.5", "config.json", "extra_arg"]
        # Function should still parse the first four arguments correctly
        expected = (10, 1.5, "config.json")
        result = validate_arguments(args)
        self.assertEqual(result, expected)

# Test class for DirectoryManager
class TestDirectoryManager(unittest.TestCase):
    def test_create_directory_if_not_exists(self):
        # Setup: Define a test directory
        test_directory = "test_directory"

        # Ensure the directory does not already exist
        if os.path.exists(test_directory):
            os.rmdir(test_directory)  # Remove if it already exists for a clean test

        # Test the function
        create_directory_if_not_exists(test_directory)

        # Assert the directory now exists
        self.assertTrue(os.path.exists(test_directory))

        # Cleanup: Remove the test directory after test
        os.rmdir(test_directory)

# Test class for save_data_to_csv
class TestSaveDataToCSV:
    """
    Test suite for the save_data_to_csv function.
    """

    @pytest.fixture
    def sample_dataframe(self):
        """
        Fixture to provide a sample DataFrame for testing.
        """
        return pd.DataFrame({
            "Column1": [1, 2, 3],
            "Column2": ["A", "B", "C"]
        })

    def test_save_data_to_csv(self, sample_dataframe, tmp_path):
        """
        Test that save_data_to_csv correctly writes a DataFrame to a CSV file.
        """
        # Create a temporary file path
        temp_file = tmp_path / "test_output.csv"
        
        # Call the function with the sample DataFrame and temporary file path
        save_data_to_csv(sample_dataframe, temp_file)
        
        # Read the file back and verify its contents
        loaded_df = pd.read_csv(temp_file)
        pd.testing.assert_frame_equal(loaded_df, sample_dataframe)

    def test_save_data_to_csv_with_mock(self, sample_dataframe):
        """
        Test that save_data_to_csv calls the correct DataFrame method using a mock.
        """
        with patch("pandas.DataFrame.to_csv") as mock_to_csv:
            # Call the function
            save_data_to_csv(sample_dataframe, "mock_path.csv")
            
            # Verify that to_csv was called once with the correct arguments
            mock_to_csv.assert_called_once_with("mock_path.csv", index=False)