import numpy as np
import pandas as pd


def generate_dummy_data_wo_interaction(
    n_cells=3000,  # cells of major cell types per individual
    sd_celltypes=0.1,  # standard deviation of number of cells
    n_major_cell_types=7,
    n_minor_cell_types=3,
    relative_abundance=0.1,  # ratio between major and rare
    n_major_diff_celltypes=1,
    n_minor_diff_celltypes=1,
    n_individuals=30,  # total individuals
    n_batches=4,
    prop_sex=0.5,
    prop_disease=0.5,
    fc_interact=0.1,  # additional proportion of interacted cells
    seed=1234,
):
    np.random.seed(seed)
    n_cell_types = n_major_cell_types + n_minor_cell_types

    # Generate unique subject IDs
    subject_id = [f"SUB_{i+1}" for i in range(n_individuals)]

    # Generate sex
    sex = np.random.choice([1, 0], size=n_individuals, p=[prop_sex, 1 - prop_sex])

    # Generate disease status
    disease = np.random.choice(
        [1, 0], size=n_individuals, p=[prop_disease, 1 - prop_disease]
    )

    # Generate age and BMI
    age = np.random.randint(18, 61, size=n_individuals)
    bmi = np.random.randint(15, 36, size=n_individuals)

    # Generate batch
    batch = np.repeat(np.arange(1, n_batches + 1), np.ceil(n_individuals / n_batches))[
        :n_individuals
    ]

    # Create dummy data DataFrame
    dummy_data = pd.DataFrame(
        {
            "subject_id": subject_id,
            "sex": pd.Categorical(sex),
            "disease": pd.Categorical(disease),
            "age": age,
            "batch": pd.Categorical(batch),
            "bmi": bmi,
        }
    )

    # Major and rare cell type counts
    major_cell_types = n_major_cell_types
    rare_cell_types = n_minor_cell_types

    celltype_df = pd.DataFrame(columns=["cell_type", "subject_id"])

    # Generate baseline of cell type data
    for id in dummy_data["subject_id"]:
        idx = dummy_data.index[dummy_data["subject_id"] == id][0]

        np.random.seed(seed * 5 * idx * 10)
        major_cell_counts = np.round(
            np.random.uniform(
                n_cells - n_cells * sd_celltypes,
                n_cells + n_cells * sd_celltypes,
                major_cell_types,
            )
        ).astype(int)

        np.random.seed(seed * 6 * idx * 10)
        rare_cell_counts = np.round(
            np.random.uniform(
                n_cells * relative_abundance
                - n_cells * relative_abundance * sd_celltypes,
                n_cells * relative_abundance
                + n_cells * relative_abundance * sd_celltypes,
                rare_cell_types,
            )
        ).astype(int)

        cell_counts = np.concatenate([major_cell_counts, rare_cell_counts])

        for i in range(n_cell_types):
            n = cell_counts[i]
            celltype_df = pd.concat(
                [
                    celltype_df,
                    pd.DataFrame(
                        {
                            "cell_type": [chr(65 + i)] * n,  # 'A' is 65 in ASCII
                            "subject_id": [id] * n,
                        }
                    ),
                ],
                ignore_index=True,
            )

    # Differential cell types
    diff_clusters = list(range(1, n_major_diff_celltypes + 1)) + list(
        range(n_cell_types - n_minor_diff_celltypes + 1, n_cell_types + 1)
    )
    diff_cell_types = [chr(65 + i) for i in diff_clusters]
    prop_control_types = 0.1

    for cell_type in [chr(65 + i) for i in range(n_cell_types)]:
        abundance = (
            celltype_df[celltype_df["cell_type"] == cell_type]
            .merge(dummy_data, on="subject_id")
            .groupby("subject_id")
            .agg(pro=("cell_type", "size"), count=("cell_type", "size"))
            .reset_index()
        )

        abundance["pro"] /= n_cells

        if cell_type in diff_cell_types:
            diff = np.ceil(
                n_cells * (abundance["pro"].median() * (1 + fc_interact))
                - n_cells * abundance["pro"].median()
            ).astype(int)
            len_disease = len(dummy_data[dummy_data["disease"] == 1])
            len_diff = len_disease * diff

            celltype_df = pd.concat(
                [
                    celltype_df,
                    pd.DataFrame(
                        {
                            "cell_type": [cell_type] * len_diff,
                            "subject_id": np.random.choice(
                                dummy_data[dummy_data["disease"] == 1]["subject_id"],
                                len_diff,
                            ),
                        }
                    ),
                ],
                ignore_index=True,
            )
        else:
            diff_control = np.ceil(
                n_cells * (abundance["pro"].median() * (1 + fc_interact))
                + n_cells * abundance["pro"].median()
            ).astype(int)
            len_control = len(dummy_data[dummy_data["disease"] == 0]) * diff_control

            celltype_df = pd.concat(
                [
                    celltype_df,
                    pd.DataFrame(
                        {
                            "cell_type": [cell_type] * len_control,
                            "subject_id": np.random.choice(
                                dummy_data[dummy_data["disease"] == 0]["subject_id"],
                                len_control,
                            ),
                        }
                    ),
                ],
                ignore_index=True,
            )

    dummy_data = dummy_data.merge(celltype_df, on="subject_id")

    # Shuffle rows
    np.random.seed(seed * 7)
    dummy_data = dummy_data.sample(frac=1).reset_index(drop=True)

    return dummy_data, diff_cell_types
