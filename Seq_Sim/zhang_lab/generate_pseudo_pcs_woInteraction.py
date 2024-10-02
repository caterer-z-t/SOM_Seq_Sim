import numpy as np
import pandas as pd


def generate_pseudo_pcs_woInteraction(
    data,
    n_pcs=20,
    cluster_pcs=range(1, 21),
    disease_pcs=0,
    sex_pcs=0,
    age_pcs=0,
    bmi_pcs=0,
    batch_pcs=0,
    scale_factor=2,
    cluster_ratio=0.25,
    disease_ratio=0,
    sex_ratio=0,
    age_ratio=0,
    bmi_ratio=0,
    batch_ratio=0,
    cluster_col="cell_type",
    disease_col="disease",
    sex_col="sex",
    age_col="age",
    bmi_col="bmi",
    batch_col="batch",
    seed=1234,
):

    np.random.seed(seed)

    # Sample the PCs
    disease_pcs = np.random.permutation(disease_pcs) if disease_pcs else []
    np.random.seed(seed * 2)
    sex_pcs = np.random.permutation(sex_pcs) if sex_pcs else []
    np.random.seed(seed * 3)
    age_pcs = np.random.permutation(age_pcs) if age_pcs else []
    np.random.seed(seed * 4)
    bmi_pcs = np.random.permutation(bmi_pcs) if bmi_pcs else []
    np.random.seed(seed * 5)
    batch_pcs = np.random.permutation(batch_pcs) if batch_pcs else []

    n_cells = data.shape[0]
    cell_clusters = pd.factorize(data[cluster_col])[0] + 1  # Start from 1 instead of 0

    pseudo_pcs = []

    for x in range(1, n_pcs + 1):
        # Shuffle clusters
        shuffled_clusters = np.random.permutation(np.unique(cell_clusters)) * 10
        for i in range(1, len(np.unique(cell_clusters)) + 1):
            cell_clusters[cell_clusters == i] = int(shuffled_clusters[i - 1])

        cell_sex = pd.factorize(data[sex_col])[0] + 1
        cell_age = pd.factorize(data[age_col])[0] + 1
        cell_bmi = pd.factorize(data[bmi_col])[0] + 1
        cell_batch = pd.factorize(data[batch_col])[0] + 1
        cell_diseases = pd.factorize(data[disease_col])[0] + 1

        variance = 1 / (x * scale_factor)

        np.random.seed(seed * x)
        pc = np.random.normal(0, np.sqrt(variance), n_cells)
        np.random.seed(seed * x)
        pc_cluster = np.random.normal(np.std(cell_clusters), np.sqrt(variance), n_cells)

        np.random.seed(seed * x)
        pc_disease = np.random.normal(np.std(cell_diseases), np.sqrt(variance), n_cells)
        np.random.seed(seed * x * 2)
        pc_sex = np.random.normal(np.std(cell_sex), np.sqrt(variance), n_cells)
        np.random.seed(seed * x * 3)
        pc_age = np.random.normal(np.std(cell_age), np.sqrt(variance), n_cells)
        np.random.seed(seed * x * 4)
        pc_bmi = np.random.normal(np.std(cell_bmi), np.sqrt(variance), n_cells)
        np.random.seed(seed * x * 5)
        pc_batch = np.random.normal(np.std(cell_batch), np.sqrt(variance), n_cells)

        # Calculate the ratios for each attribute
        cluster_ratio_tmp = cluster_ratio / (x if x in cluster_pcs else 1)
        disease_ratio_tmp = disease_ratio / (x if x in disease_pcs else 1)
        sex_ratio_tmp = sex_ratio / (x if x in sex_pcs else 1)
        age_ratio_tmp = age_ratio / (x if x in age_pcs else 1)
        bmi_ratio_tmp = bmi_ratio / (x if x in bmi_pcs else 1)
        batch_ratio_tmp = batch_ratio / (x if x in batch_pcs else 1)

        pseudo_pc = (
            pc
            * (
                1
                - cluster_ratio_tmp
                - disease_ratio_tmp
                - sex_ratio_tmp
                - age_ratio_tmp
                - bmi_ratio_tmp
                - batch_ratio_tmp
            )
            + pc_cluster * cluster_ratio_tmp
            + pc_disease * disease_ratio_tmp
            + pc_sex * sex_ratio_tmp
            + pc_age * age_ratio_tmp
            + pc_bmi * bmi_ratio_tmp
            + pc_batch * batch_ratio_tmp
        )

        pseudo_pcs.append(pseudo_pc)

    return np.array(pseudo_pcs).T  # Transpose to get shape (n_cells, n_pcs)
