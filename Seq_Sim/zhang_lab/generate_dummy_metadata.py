import numpy as np
import pandas as pd


def generate_dummy_metadata(
    n_cells=3000,  # cells of major cell types per individual
    sd_celltypes=0.1,  # standard deviation of number of cells
    n_major_cell_types=7,
    n_minor_cell_types=3,
    relative_abundance=0.1,  # ratio between major and rare
    n_major_diff_celltypes=1,
    n_minor_diff_celltypes=1,
    n_individuals=30,  # total individuals
    n_batchs=4,
    prop_sex=0.5,
    prop_disease=0.5,
    fc_increase=0.1,  # additional FC of specified cell types which are from people with case groups compared to control groups
    seed=1234,
):

    n_cell_types = n_major_cell_types + n_minor_cell_types

    # Generate unique subject IDs
    np.random.seed(seed)
    subject_id = [f"SUB_{i+1}" for i in range(n_individuals)]

    sex = np.random.choice([1, 0], size=n_individuals, p=[prop_sex, 1 - prop_sex])

    np.random.seed(seed * 2)
    disease = np.random.choice(
        [1, 0], size=n_individuals, p=[prop_disease, 1 - prop_disease]
    )

    np.random.seed(seed * 3)
    age = np.random.randint(18, 61, size=n_individuals)

    np.random.seed(seed * 4)
    bmi = np.random.randint(15, 36, size=n_individuals)

    batch = np.tile(np.arange(1, n_batchs + 1), int(np.ceil(n_individuals / n_batchs)))[
        :n_individuals
    ]

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

    major_cell_types = n_major_cell_types
    rare_cell_types = n_minor_cell_types

    celltype_df = pd.DataFrame(columns=["cell_type", "subject_id"])

    # Generate baseline of cell type data.frame
    for id in dummy_data["subject_id"]:
        np.random.seed(seed * 5 * np.where(dummy_data["subject_id"] == id)[0][0] * 10)
        major_cell_counts = np.round(
            np.random.uniform(
                n_cells - n_cells * sd_celltypes,
                n_cells + n_cells * sd_celltypes,
                major_cell_types,
            )
        ).astype(int)

        np.random.seed(seed * 6 * np.where(dummy_data["subject_id"] == id)[0][0] * 10)
        rare_cell_counts = np.round(
            np.random.uniform(
                n_cells * relative_abundance
                - n_cells * relative_abundance * sd_celltypes,
                n_cells * relative_abundance
                + n_cells * relative_abundance * sd_celltypes,
                rare_cell_types,
            )
        ).astype(int)

        cell_counts = np.concatenate((major_cell_counts, rare_cell_counts))

        for i in range(n_cell_types):
            n = cell_counts[i]
            celltype_df = pd.concat(
                [
                    celltype_df,
                    pd.DataFrame(
                        {"cell_type": [chr(65 + i)] * n, "subject_id": [id] * n}
                    ),
                ]
            )

    # Identifying different cell types
    diff_clusters = list(range(1, n_major_diff_celltypes + 1)) + list(
        range(n_cell_types - n_minor_diff_celltypes + 1, n_cell_types + 1)
    )
    diff_cell_types = [chr(65 + cluster - 1) for cluster in diff_clusters]

    for i in range(n_cell_types):
        cell_type_letter = chr(65 + i)
        if cell_type_letter in diff_cell_types:
            abundance = celltype_df[celltype_df["cell_type"] == cell_type_letter].merge(
                dummy_data, on="subject_id"
            )
            pro = abundance["cell_type"].count() / n_cells
            diff = int(np.ceil(n_cells * (pro * (1 + fc_increase)) - n_cells * pro))
            len_diff = (
                abundance[dummy_data["disease"] == 1]["subject_id"].nunique() * diff
            )

            celltype_df = pd.concat(
                [
                    pd.DataFrame(
                        {
                            "cell_type": [cell_type_letter] * len_diff,
                            "subject_id": np.repeat(
                                abundance[dummy_data["disease"] == 1][
                                    "subject_id"
                                ].unique(),
                                diff,
                            ),
                        }
                    ),
                    celltype_df,
                ]
            )

    dummy_data = dummy_data.merge(celltype_df, on="subject_id", how="left")

    # Shuffle rows
    np.random.seed(seed * 7)
    dummy_data = dummy_data.sample(frac=1).reset_index(drop=True)

    return dummy_data, diff_cell_types
