import numpy as np
import pandas as pd
from joblib import Parallel, delayed


def generate_pseudo_features(
    data,
    n_features=20,
    cluster_features=list(range(1, 21)),
    disease_features=0,
    sex_features=0,
    age_features=0,
    bmi_features=0,
    batch_features=0,
    individual_features=0,
    cluster_ratio=0.25,
    disease_ratio=0,
    sex_ratio=0,
    age_ratio=0,
    bmi_ratio=0,
    batch_ratio=0,
    individual_ratio=0.1,
    ratio_variance=0.5,
    cluster_col="cell_type",
    disease_col="disease",
    sex_col="sex",
    age_col="age",
    bmi_col="bmi",
    batch_col="batch",
    individual_col="individual",
    seed=1234,
    n_thread=1,
    verbose=True,
):

    np.random.seed(seed)

    # Sampling features (assuming features are given in an iterable)
    disease_features = np.random.permutation(disease_features)
    sex_features = np.random.permutation(sex_features)
    age_features = np.random.permutation(age_features)
    bmi_features = np.random.permutation(bmi_features)
    batch_features = np.random.permutation(batch_features)
    individual_features = np.random.permutation(individual_features)

    n_cells = data.shape[0]
    n_clusters = data[cluster_col].nunique()

    cell_sex = data[sex_col].astype("category").cat.codes
    cell_age = data[age_col].astype("category").cat.codes
    cell_bmi = data[bmi_col].astype("category").cat.codes
    cell_batch = data[batch_col].astype("category").cat.codes
    cell_diseases = data[disease_col].astype("category").cat.codes
    cell_individual = data[individual_col].astype("category").cat.codes

    def generate_feature(x):
        var_all = []

        # Generate dummy features reflecting cell clusters
        cell_clusters = data[cluster_col].astype("category").cat.codes
        for j in range(1, 3):  # 2 iterations
            cell_clusters_tmp = cell_clusters.copy()
            np.random.seed(x * j)
            cluster_mean = (
                np.random.choice(np.unique(cell_clusters_tmp), size=n_clusters) * 10
            )
            for i in range(n_clusters):
                cell_clusters_tmp[cell_clusters == i] = int(cluster_mean[i])
            if verbose:
                print(f"length of cell_clusters_tmp: {len(cell_clusters_tmp)}")
        variance = 1 / cluster_mean  # cell type specific variance
        np.random.seed(seed * x)
        pc_cluster = np.random.normal(
            loc=np.mean(cell_clusters), scale=np.sqrt(variance), size=n_cells
        )
        var_all.append(variance)

        # Similar process for disease, sex, age, bmi, batch, and individual
        for covariate, covariate_col in zip(
            [cell_diseases, cell_sex, cell_age, cell_bmi, cell_batch, cell_individual],
            [disease_col, sex_col, age_col, bmi_col, batch_col, individual_col],
        ):
            for j in range(1, 3):
                covariates_tmp = covariate.copy()
                np.random.seed(
                    x
                    * j
                    * (
                        covariate_col == disease_col
                        and 2
                        or covariate_col == sex_col
                        and 3
                        or covariate_col == age_col
                        and 4
                        or covariate_col == bmi_col
                        and 5
                        or covariate_col == batch_col
                        and 6
                        or 7
                    )
                )
                cluster_mean = (
                    np.random.choice(
                        np.unique(covariates_tmp), size=len(np.unique(covariates_tmp))
                    )
                    * 10
                )
                for i in range(len(np.unique(covariate))):
                    covariates_tmp[covariate == i] = int(cluster_mean[i])
            variance = 1 / cluster_mean  # cell type specific variance
            np.random.seed(
                seed
                * x
                * (
                    covariate_col == disease_col
                    and 2
                    or covariate_col == sex_col
                    and 3
                    or covariate_col == age_col
                    and 4
                    or covariate_col == bmi_col
                    and 5
                    or covariate_col == batch_col
                    and 6
                    or 7
                )
            )
            pc_covariate = np.random.normal(
                loc=np.mean(covariates_tmp), scale=np.sqrt(variance), size=n_cells
            )
            var_all.append(variance)

        return var_all

    # Parallel processing
    results = Parallel(n_jobs=n_thread)(
        delayed(generate_feature)(x) for x in range(1, n_features + 1)
    )

    # Combine results (depending on how you want to store these features)
    combined_results = pd.DataFrame(results)

    return combined_results
