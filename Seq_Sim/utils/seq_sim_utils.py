import numpy as np
import pandas as pd
from typing import List, Tuple, Dict, Any
import yaml
import os
import argparse


# Constant for seed to ensure reproducibility
DEFAULT_SEED = 1234

# Function to generate dummy data with interactions
def generate_feature(idx, data, seed, cluster_ratio, ratio_variance, 
                     cluster_col, disease_col, individual_col):
    """Process a single feature for pseudo-feature generation.

    Args:
        idx (int): Feature index.
        data (pd.DataFrame): Dummy data.
        seed (int): Seed for reproducibility.
        cluster_ratio (float): Cluster ratio.
        ratio_variance (float): Ratio variance.
        cluster_col (str): Cluster column name.
        disease_col (str): Disease column name.
        individual_col (str): Individual column name.

    Returns:
        np.ndarray: Generated feature.
    """
    return process_feature(idx, data, seed, cluster_ratio, ratio_variance, 
                           cluster_col, disease_col, individual_col)

# Function to generate dummy data from pseudo-features generation
def generate_all_features(n_features, data, seed, cluster_ratio, 
                          ratio_variance, cluster_col, disease_col, 
                          individual_col):
    """Generate all features for pseudo-feature generation.

    Args:
        n_features (int): Number of features to generate.
        data (pd.DataFrame): Dummy data.
        seed (int): Seed for reproducibility.
        cluster_ratio (float): Cluster ratio.
        ratio_variance (float): Ratio variance.
        cluster_col (str): Cluster column name.
        disease_col (str): Disease column name.
        individual_col (str): Individual column name.

    Returns:
        List of generated features.
    """
    # Initialize list to store features
    features = []

    # Generate features
    for idx in range(1, n_features + 1):
        feature = generate_feature(idx, data, seed, cluster_ratio, 
                                   ratio_variance, cluster_col, 
                                   disease_col, individual_col)
        features.append(feature)
    return features

# Function to combine features into a DataFrame
def combine_features(features, n_features):
    """
    Combine the features into a DataFrame.

    Args:
        features (List): List of features.
        n_features (int): Number of features.

    Returns:
        pd.DataFrame: Pseudo-features
    """
    # Create a DataFrame from the features
    features_df = pd.DataFrame(features).T

    # Rename the columns
    features_df.columns = [f"feature_{i}" for i in range(1, n_features + 1)]
    return features_df

# Function to generate pseudo-features for simulation
def generate_pseudo_features(data,
                             n_features=20,
                             cluster_ratio=0.25,
                             ratio_variance=0.5,
                             cluster_col="cell_type",
                             disease_col="disease",
                             individual_col="batch",
                             seed=1234):
    """
    Generate pseudo-features for simulation.

    Args:
        data (pd.DataFrame): Dummy data.
        n_features (int, optional): Number of features to generate. Defaults to 20.
        cluster_ratio (float, optional): Cluster ratio. Defaults to 0.25.
        ratio_variance (float, optional): Ratio variance. Defaults to 0.5.
        cluster_col (str, optional): Cluster column. Defaults to "cell_type".
        disease_col (str, optional): Disease column. Defaults to "disease".
        individual_col (str, optional): Individual column. Defaults to "batch".
        seed (int, optional): Seed. Defaults to 1234.

    Returns:
        pd.DataFrame: Pseudo-features.
    """
    # Step 1: Generate all features
    features = generate_all_features(n_features, data, seed, cluster_ratio, 
                                     ratio_variance, cluster_col, 
                                     disease_col, individual_col)

    # Step 2: Combine the features into a DataFrame
    features_df = combine_features(features, n_features)
    
    return features_df

# Function to set random seed for reproducibility
def set_random_seed(x, seed):
    """Set random seed for reproducibility.

    Args:
        x (int): Feature index.
        seed (int): Seed for reproducibility.
    """
    np.random.seed(x * seed)

# Function to encode categorical columns for pseudo-feature generation
def encode_categorical_columns(data, cluster_col, disease_col, individual_col):
    """
    Encode categorical columns for pseudo-feature generation.

    Args:
        data (pd.DataFrame): Dummy data.
        cluster_col (str): Cluster column name.
        disease_col (str): Disease column name.
        individual_col (str): Individual column name.

    Returns:
        Tuple of encoded cell diseases, cell individuals, and cell clusters.
    """

    # Encode categorical columns
    cell_diseases = pd.factorize(data[disease_col])[0]
    cell_individual = pd.factorize(data[individual_col])[0]
    cell_clusters = pd.factorize(data[cluster_col])[0]
    
    return cell_diseases, cell_individual, cell_clusters

# Function to generate features based on cluster data
def generate_variance(data):
    """
    Generate variance for the given data.

    Args:
        data (np.ndarray): The data based on which variance is to be generated.

    Returns:
        float: Mean of the generated variance.
    """
    return np.mean(1 / (data + 1e-8))

# Function to generate cluster features
def generate_cluster_features(cell_clusters, unique_clusters):
    """ 
    Generate features based on cluster data.

    Args:
        cell_clusters (np.ndarray): Encoded cell clusters.
        unique_clusters (int): Number of unique clusters.

    Returns:
        Tuple of cluster features and variance.
    """
    # Generate cluster features
    cluster_mean = np.random.choice(unique_clusters, unique_clusters) * 10
    cell_clusters_tmp = np.array([cluster_mean[cell_clusters[i]] for i in range(len(cell_clusters))])
    
    # Generate variance
    variance = generate_variance(cell_clusters_tmp)
    
    # Generate final cluster features with noise
    pc_cluster = np.random.normal(loc=cell_clusters_tmp, scale=np.sqrt(variance), size=len(cell_clusters))
    
    # Return cluster features and variance
    return pc_cluster, variance

# Function to generate disease features
def generate_disease_features(cell_diseases):
    """ 
    Generate features based on disease data.

    Args:
        cell_diseases (np.ndarray): Encoded cell diseases.

    Returns:
        float: Variance of disease features.
    """
    # Generate variance using the common function
    return generate_variance(cell_diseases)

# Function to generate individual features
def generate_individual_features(cell_individual):
    """Generate features based on individual data.

    Args:
        cell_individual (np.ndarray): Encoded cell individuals.

    Returns:
        float: Variance of individual features.
    """
    # Generate variance using the common function
    return generate_variance(cell_individual)

# Function to apply noise to the generated features
def apply_noise(pc_cluster, var_all, n_cells, cluster_ratio, ratio_variance):
    """Apply noise to the generated features.

    Args:
        pc_cluster (int): Cluster features.
        var_all (List): List of variances.
        n_cells (int): Number of cells.
        cluster_ratio (float): Cluster ratio.
        ratio_variance (float): Ratio variance.

    Returns:
        np.ndarray: Final feature with noise.
    """
    # Apply noise
    final_variance = np.random.choice(var_all, size=n_cells)
    noise_pc = np.random.normal(loc=0, scale=np.sqrt(final_variance), size=n_cells)
    cluster_ratio_tmp = cluster_ratio + np.random.uniform(-cluster_ratio * ratio_variance, cluster_ratio * ratio_variance)
    
    # Combine features with noise
    final_feature = pc_cluster * cluster_ratio_tmp + noise_pc
    return final_feature

# Function to process a single feature for pseudo-feature generation
def process_feature(x, data, seed, cluster_ratio, ratio_variance, cluster_col, disease_col, individual_col):
    """
    Process a single feature for pseudo-feature generation.

    Args:
        x (int): Feature index.
        data (pd.DataFrame): Dummy data.
        seed (int): Seed for reproducibility.
        cluster_ratio (float): Cluster ratio.
        ratio_variance (float): Ratio variance.
        cluster_col (str): Cluster column name.
        disease_col (str): Disease column name.
        individual_col (str): Individual column name.

    Returns:
        np.ndarray: Generated feature.
    """
    # Step 1: Set random seed for reproducibility
    set_random_seed(x, seed)
    
    n_cells = len(data)
    unique_clusters = data[cluster_col].nunique()

    # Step 2: Encode categorical columns
    cell_diseases, cell_individual, cell_clusters = encode_categorical_columns(data, cluster_col, disease_col, individual_col)
    
    # Step 3: Generate cluster features
    pc_cluster, cluster_variance = generate_cluster_features(cell_clusters, unique_clusters)

    # Step 4: Generate disease features
    disease_variance = generate_disease_features(cell_diseases)

    # Step 5: Generate individual features
    individual_variance = generate_individual_features(cell_individual)

    # Step 6: Combine features with noise
    var_all = [cluster_variance, disease_variance, individual_variance]
    final_feature = apply_noise(pc_cluster, var_all, n_cells, cluster_ratio, ratio_variance)

    return final_feature

# Function to generate dummy data with interactions
def generate_dummy_data_wo_interaction(
    n_cells: int = 3000,
    sd_celltypes: float = 0.1,
    n_major_cell_types: int = 7,
    n_minor_cell_types: int = 3,
    relative_abundance: float = 0.1,
    n_major_diff_celltypes: int = 1,
    n_minor_diff_celltypes: int = 1,
    n_individuals: int = 30,
    n_batches: int = 4,
    prop_sex: float = 0.5,
    prop_disease: float = 0.5,
    fc_interact: float = 0.1,
    seed: int = DEFAULT_SEED
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Generate dummy data for simulation without interactions.

    Args:
        n_cells: Number of cells of major cell types per individual.
        sd_celltypes: Standard deviation of number of cells.
        n_major_cell_types: Number of major cell types.
        n_minor_cell_types: Number of minor cell types.
        relative_abundance: Ratio between major and rare cell types.
        n_major_diff_celltypes: Number of differentially expressed major cell types.
        n_minor_diff_celltypes: Number of differentially expressed minor cell types.
        n_individuals: Total number of individuals.
        n_batches: Number of batches.
        prop_sex: Proportion of sex.
        prop_disease: Proportion of disease.
        fc_interact: Proportion of interacted cells.
        seed: Random seed for reproducibility.

    Returns:
        A tuple containing:
            - DataFrame: The generated dummy data.
            - List of differential cell types.
    """
    np.random.seed(seed)
    subject_ids = [f"SUB_{i}" for i in range(1, n_individuals + 1)]
    
    sex = np.random.choice([0, 1], size=n_individuals, p=[1 - prop_sex, prop_sex])
    disease = np.random.choice([0, 1], size=n_individuals, p=[1 - prop_disease, prop_disease])
    age = np.random.randint(18, 61, size=n_individuals)
    bmi = np.random.randint(15, 36, size=n_individuals)
    batch = np.tile(range(1, n_batches + 1), int(np.ceil(n_individuals / n_batches)))[:n_individuals]

    dummy_data = pd.DataFrame({
        "subject_id": subject_ids,
        "sex": pd.Categorical(sex, categories=[0, 1]),
        "disease": pd.Categorical(disease, categories=[0, 1]),
        "age": age,
        "batch": pd.Categorical(batch, categories=range(1, n_batches + 1)),
        "bmi": bmi
    })
    
    # Generate baseline cell counts
    major_cell_counts, rare_cell_counts = generate_cell_counts(
        subject_ids, n_cells, sd_celltypes, n_major_cell_types, n_minor_cell_types, relative_abundance
    )
    
    # Create cell type dataframe
    celltype_df = create_celltype_dataframe(
        subject_ids, major_cell_counts, rare_cell_counts, n_major_cell_types, n_minor_cell_types
    )
    
    # Add differential cell types
    diff_cell_types = add_differential_cell_types(
        celltype_df, dummy_data, n_cells, n_major_diff_celltypes, n_minor_diff_celltypes, fc_interact
    )
    
    # Merge and shuffle the data
    dummy_data = dummy_data.merge(celltype_df, on="subject_id", how="left")
    dummy_data = dummy_data.sample(frac=1, random_state=seed).reset_index(drop=True)

    return dummy_data, diff_cell_types

# Function to generate major cell counts
def generate_major_cell_counts(n_cells: int, sd_celltypes: float, n_major_cell_types: int) -> List[int]:
    """
    Generate cell counts for major cell types based on a uniform distribution.

    Args:
        n_cells: Total number of cells per subject.
        sd_celltypes: Standard deviation of cell counts.
        n_major_cell_types: Number of major cell types.

    Returns:
        List of major cell counts.
    """
    return np.random.uniform(n_cells * (1 - sd_celltypes), n_cells * (1 + sd_celltypes), n_major_cell_types).round()

# Function to generate rare cell counts
def generate_rare_cell_counts(n_cells: int, sd_celltypes: float, relative_abundance: float, n_minor_cell_types: int) -> List[int]:
    """
    Generate cell counts for minor cell types based on a uniform distribution.

    Args:
        n_cells: Total number of cells per subject.
        sd_celltypes: Standard deviation of cell counts.
        relative_abundance: Ratio between major and minor cell types.
        n_minor_cell_types: Number of minor cell types.

    Returns:
        List of rare cell counts.
    """
    return np.random.uniform(n_cells * relative_abundance * (1 - sd_celltypes),
                             n_cells * relative_abundance * (1 + sd_celltypes), n_minor_cell_types).round()

# Function to generate cell counts for a single subject
def generate_cell_counts_for_subject(
    n_cells: int,
    sd_celltypes: float,
    n_major_cell_types: int,
    n_minor_cell_types: int,
    relative_abundance: float
) -> Tuple[List[int], List[int]]:
    """
    Generate major and rare cell counts for a single subject.

    Args:
        n_cells: Total number of cells per subject.
        sd_celltypes: Standard deviation of cell counts.
        n_major_cell_types: Number of major cell types.
        n_minor_cell_types: Number of minor cell types.
        relative_abundance: Ratio between major and minor cell types.

    Returns:
        Tuple of major cell counts and rare cell counts.
    """
    # Generate major and rare cell counts
    major_cell_counts = generate_major_cell_counts(n_cells, sd_celltypes, n_major_cell_types)
    rare_cell_counts = generate_rare_cell_counts(n_cells, sd_celltypes, relative_abundance, n_minor_cell_types)
    
    return major_cell_counts, rare_cell_counts

# Function to generate baseline cell counts
def generate_cell_counts(
    subject_ids: List[str],
    n_cells: int,
    sd_celltypes: float,
    n_major_cell_types: int,
    n_minor_cell_types: int,
    relative_abundance: float
) -> Tuple[List[List[int]], List[List[int]]]:
    """
    Generate baseline cell counts for major and minor cell types for all subjects.

    Args:
        subject_ids: List of subject IDs.
        n_cells: Total number of cells per subject.
        sd_celltypes: Standard deviation of cell counts.
        n_major_cell_types: Number of major cell types.
        n_minor_cell_types: Number of minor cell types.
        relative_abundance: Ratio between major and minor cell types.

    Returns:
        A tuple of lists: major cell counts and rare cell counts for each subject.
    """
    # Initialize lists to store cell counts
    major_cell_counts_all_subjects = []
    rare_cell_counts_all_subjects = []
    
    # Generate cell counts for each subject
    for _ in subject_ids:
        major_cell_counts, rare_cell_counts = generate_cell_counts_for_subject(
            n_cells, sd_celltypes, n_major_cell_types, n_minor_cell_types, relative_abundance
        )
        major_cell_counts_all_subjects.append(major_cell_counts)
        rare_cell_counts_all_subjects.append(rare_cell_counts)

    # Return the cell counts for all subjects
    return major_cell_counts_all_subjects, rare_cell_counts_all_subjects

# Function to create new rows for a specific cell type
def create_new_rows_for_cell_type(subj_id: str, cell_type: str, count: int) -> pd.DataFrame:
    """
    Create new rows for a specific cell type for a given subject.

    Args:
        subj_id: The subject ID.
        cell_type: The type of cell.
        count: The number of cells of that type.

    Returns:
        DataFrame with new rows for the given cell type.
    """
    return pd.DataFrame({
        "cell_type": [cell_type] * int(count),
        "subject_id": [subj_id] * int(count)
    })

# Function to create the dataframe for cell types
def create_celltype_dataframe(
    subject_ids: List[str],
    major_cell_counts: List[List[int]],
    rare_cell_counts: List[List[int]],
    n_major_cell_types: int,
    n_minor_cell_types: int
) -> pd.DataFrame:
    """
    Create the dataframe for cell types based on the given counts.

    Args:
        subject_ids: List of subject IDs.
        major_cell_counts: Major cell counts for each subject.
        rare_cell_counts: Rare cell counts for each subject.
        n_major_cell_types: Number of major cell types.
        n_minor_cell_types: Number of minor cell types.

    Returns:
        A DataFrame containing the cell types for each subject.
    """
    # Step 1: Generate the list of cell types
    cell_types = generate_cell_types(n_major_cell_types, n_minor_cell_types)

    # Step 2: Create an empty DataFrame to store the cell type data
    celltype_df = pd.DataFrame(columns=["cell_type", "subject_id"])

    # Step 3: Loop through subjects and their respective cell counts
    for i, subj_id in enumerate(subject_ids):
        counts = list(major_cell_counts[i]) + list(rare_cell_counts[i])
        for j, count in enumerate(counts):
            # Step 4: Create new rows for each cell type and subject
            new_rows = create_new_rows_for_cell_type(subj_id, cell_types[j], count)
            # Step 5: Append the new rows to the DataFrame
            celltype_df = pd.concat([celltype_df, new_rows], ignore_index=True)

    return celltype_df

# Function to identify differentially expressed clusters
def identify_diff_clusters(n_major_diff_celltypes: int, n_minor_diff_celltypes: int, n_cells: int) -> List[int]:
    """
    Identify the indices of the differentially expressed clusters.

    Args:
        n_major_diff_celltypes: Number of differentially expressed major cell types.
        n_minor_diff_celltypes: Number of differentially expressed minor cell types.
        n_cells: Total number of cells per individual.

    Returns:
        List of indices representing the differentially expressed clusters.
    """
    return list(range(n_major_diff_celltypes)) + list(range(n_cells - n_minor_diff_celltypes, n_cells))

# Helper function to generate cell types
def generate_cell_types_from_range(start_index: int, num_types: int) -> List[str]:
    """
    Generate a list of cell types based on a starting index and number of types.

    Args:
        start_index: Starting index for the cell types.
        num_types: Number of cell types to generate.

    Returns:
        List of cell types.
    """
    return [chr(65 + start_index + i) for i in range(num_types)]

# Function to generate differential cell types (reuses generate_cell_types_from_range)
def generate_diff_cell_types(diff_clusters: List[int]) -> List[str]:
    """
    Generate the list of differential cell types based on the clusters.

    Args:
        diff_clusters: List of indices for the differentially expressed clusters.

    Returns:
        List of differential cell types.
    """
    return generate_cell_types_from_range(0, len(diff_clusters))

# Function to generate cell types (reuses generate_cell_types_from_range)
def generate_cell_types(n_major_cell_types: int, n_minor_cell_types: int) -> List[str]:
    """
    Generate a list of cell types based on the number of major and minor cell types.

    Args:
        n_major_cell_types: Number of major cell types.
        n_minor_cell_types: Number of minor cell types.

    Returns:
        List of cell types.
    """
    return generate_cell_types_from_range(0, n_major_cell_types + n_minor_cell_types)

# Function to calculate the abundance of a specific cell type
def calculate_abundance(celltype_df: pd.DataFrame, dummy_data: pd.DataFrame, cell_type: str, n_cells: int) -> np.ndarray:
    """
    Calculate the abundance of a specific cell type across individuals.

    Args:
        celltype_df: DataFrame of cell types.
        dummy_data: DataFrame of dummy data.
        cell_type: The specific cell type to calculate abundance for.
        n_cells: Total number of cells per individual.

    Returns:
        Abundance of the cell type for each subject.
    """
    abundance = (
        celltype_df[celltype_df["cell_type"] == cell_type]
        .merge(dummy_data, on="subject_id")
        .groupby("subject_id")
        .size()
        / n_cells
    )
    return abundance

# Function to calculate differential expression based on fold change
def calculate_diff_expression(abundance: np.ndarray, n_cells: int, fc_interact: float) -> int:
    """
    Calculate the differential expression based on fold change.

    Args:
        abundance: Abundance of the cell type for each subject.
        n_cells: Total number of cells per individual.
        fc_interact: Fold change for interacted cells.

    Returns:
        The number of additional rows to add to the dataframe for the cell type.
    """
    diff = np.ceil(n_cells * (np.median(abundance) * (1 + fc_interact)) - n_cells * np.median(abundance)).astype(int)
    return diff

# Function to add additional rows for the differentially expressed cell type
def add_rows_for_diff_cells(celltype_df: pd.DataFrame, dummy_data: pd.DataFrame, diff: int, cell_type: str) -> pd.DataFrame:
    """
    Add additional rows to the dataframe for the differentially expressed cell type.

    Args:
        celltype_df: DataFrame of cell types.
        dummy_data: DataFrame of dummy data.
        diff: Number of additional rows to add for the cell type.
        cell_type: The cell type for which to add rows.

    Returns:
        The updated dataframe with added rows for the differential cell type.
    """

    # Add rows for the differentially expressed cell type
    for subj_id in dummy_data.loc[dummy_data["disease"] == 1, "subject_id"]:
        additional_rows = pd.DataFrame({
            "cell_type": [cell_type] * diff,
            "subject_id": [subj_id] * diff
        })
        celltype_df = pd.concat([celltype_df, additional_rows], ignore_index=True)

    return celltype_df

# Function to add differential expression for cell types based on disease condition
def add_differential_cell_types(
    celltype_df: pd.DataFrame,
    dummy_data: pd.DataFrame,
    n_cells: int,
    n_major_diff_celltypes: int,
    n_minor_diff_celltypes: int,
    fc_interact: float
) -> List[str]:
    """
    Add differential expression for cell types based on disease condition.

    Args:
        celltype_df: DataFrame of cell types.
        dummy_data: DataFrame of dummy data.
        n_cells: Total number of cells per individual.
        n_major_diff_celltypes: Number of differentially expressed major cell types.
        n_minor_diff_celltypes: Number of differentially expressed minor cell types.
        fc_interact: Fold change for interacted cells.

    Returns:
        List of differential cell types.
    """
    # Step 1: Identify differential clusters
    diff_clusters = identify_diff_clusters(n_major_diff_celltypes, n_minor_diff_celltypes, n_cells)
    
    # Step 2: Generate differential cell types
    diff_cell_types = generate_diff_cell_types(diff_clusters)
    
    # Step 3: For each differential cell type, calculate abundance and add rows
    for cell_type in diff_cell_types:
        abundance = calculate_abundance(celltype_df, dummy_data, cell_type, n_cells)
        
        # Step 4: Calculate differential expression
        diff = calculate_diff_expression(abundance, n_cells, fc_interact)
        
        # Step 5: Add rows for the differential cell type
        celltype_df = add_rows_for_diff_cells(celltype_df, dummy_data, diff, cell_type)

    return diff_cell_types

# Function to load configuration from a YAML file
def load_config(config_file: str) -> dict:
    """Load configuration settings from a YAML file.

    Args:
        config_file: Path to the YAML configuration file.

    Returns:
        A dictionary of configuration settings.
    """
    with open(config_file, 'r') as file:
        config = yaml.safe_load(file)
    return config

# Case when function to mimic dplyr's case_when
def case_when(condition_list: List[Tuple[bool, any]], default_value: any) -> any:
    """
    Mimic dplyr's case_when functionality.

    Args:
        condition_list: List of tuples with boolean conditions and associated return values.
        default_value: The value returned if no condition matches.

    Returns:
        The value associated with the first matching condition, or the default value.
    """
    for condition, result in condition_list:
        if condition:
            return result
    return default_value


# Function to create a directory if it does not exist
def create_directory_if_not_exists(directory: str) -> None:
    """
    Creates a directory if it does not already exist.

    Args:
        directory (str): Path to the directory.
    """
    if not os.path.exists(directory):
        os.makedirs(directory)

# Function to save data to a CSV file
def save_data_to_csv(data: pd.DataFrame, file_path: str) -> None:
    """
    Saves a DataFrame to a CSV file.

    Args:
        data (pd.DataFrame): DataFrame to save.
        file_path (str): Destination file path.
    """
    data.to_csv(file_path, index=False)

# Function to generate and save features
def generate_and_save_features(num_samples: int, fold_change: float, config: Dict[str, Any]) -> None:
    """
    Generates dummy data and pseudo features, and saves them based on the configuration.

    Args:
        num_samples (int): Number of samples to generate.
        fold_change (float): Fold change for interaction effects.
        config (Dict[str, Any]): Configuration dictionary.
    """
    # Generate dummy data
    dummy_data, _ = generate_dummy_data_wo_interaction(
        n_cells=config["dummy_dataset_params"]["n_cells"],
        sd_celltypes=config["dummy_dataset_params"]["sd_celltypes"],
        n_major_cell_types=config["dummy_dataset_params"]["n_major_cell_types"],
        n_minor_cell_types=config["dummy_dataset_params"]["n_minor_cell_types"],
        relative_abundance=config["dummy_dataset_params"]["relative_abundance"],
        n_major_diff_celltypes=config["dummy_dataset_params"]["n_major_diff_celltypes"],
        n_minor_diff_celltypes=config["dummy_dataset_params"]["n_minor_diff_celltypes"],
        n_individuals=num_samples,
        n_batches=config["dummy_dataset_params"]["n_batchs"],
        prop_sex=config["dummy_dataset_params"]["prop_sex"],
        prop_disease=config["dummy_dataset_params"]["prop_disease"],
        fc_interact=fold_change,
        seed=config["dummy_dataset_params"]["seed"]
    )

    # Generate pseudo feature matrix
    pseudo_feature_matrix = generate_pseudo_features(
        data=dummy_data,
        n_features=config["dummy_dataset_params"]["n_features"],
        cluster_ratio=config["variance_attributes"]["cluster_ratio"],
        ratio_variance=config["ratio_variance"],
        cluster_col=config["column_information"]["cluster_col"],
        disease_col=config["column_information"]["disease_col"],
        individual_col=config["column_information"]["individual_col"],
        seed=config["dummy_dataset_params"]["seed"]
    )

    # Assign feature names
    pseudo_feature_matrix.columns = [f"Feature{i + 1}" for i in range(pseudo_feature_matrix.shape[1])]

    # Ensure output directory exists
    create_directory_if_not_exists(config["data_file_path"])

    # Save pseudo feature matrix if enabled
    if config["files_to_save"]["feature_matrix"]:
        pseudo_feature_file_name = os.path.join(
            config["data_file_path"],
            f"{config['file_prefix']}_pseudo_feature_num_samples_{num_samples}_fc_{fold_change}.csv"
        )
        save_data_to_csv(pseudo_feature_matrix, pseudo_feature_file_name)

    # Save latent factors if enabled
    if config["files_to_save"]["latent_factors"]:
        meta_file_name = os.path.join(
            config["data_file_path"],
            f"{config['file_prefix']}_latent_data_num_samples_{num_samples}_fc_{fold_change}.csv"
        )
        save_data_to_csv(dummy_data, meta_file_name)

def arg_parser():
    """Argument parser for the command-line arguments.

    Returns:
        argparse.ArgumentParser: Argument parser object.
    """
    parser = argparse.ArgumentParser(
        description="Generate and save features for the given number of samples and fold change."
    )
    parser.add_argument(
        "--num_samples", type=int, required=True, help="Number of samples to generate."
    )
    parser.add_argument(
        "--fold_change", type=float, required=True, help="Fold change value."
    )
    parser.add_argument(
        "--config_file", type=str, required=True, help="Path to configuration file."
    )
    return parser.parse_args()


# Function to validate and parse command-line arguments
def validate_arguments(args: list[str]) -> tuple[int, float, str]:
    """
    Validates and parses command-line arguments.

    Args:
        args (list[str]): List of command-line arguments.

    Returns:
        tuple[int, float, str]: Parsed number of samples, fold change, and config file path.

    Raises:
        ValueError: If arguments are invalid or insufficient.
    """
    if len(args) < 4:
        raise ValueError("Usage: python my_script.py <num_samples> <fold_change> <config_file>")
    try:
        num_samples = int(args[1])
        fold_change = float(args[2])
        config_file = args[3]
    except ValueError as e:
        raise ValueError("Invalid argument type. num_samples should be int and fold_change should be float.") from e
    return num_samples, fold_change, config_file
