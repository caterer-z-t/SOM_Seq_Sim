import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.preprocessing import normalize
from sklearn.linear_model import LinearRegression
from sklearn.utils.validation import check_is_fitted
from moments import kurtosis


def create_nam(
    metadata,
    pcs,
    samplem_key=None,
    graph_use="RNA_snn",
    batches=None,
    covs=None,
    nsteps=None,
    verbose=True,
    assay=None,
    key="NAMPC_",
):
    # Prepare metadata
    meta = metadata.copy()
    meta.index = range(meta.shape[0])

    m = sparse.csr_matrix(pcs.T)  # Convert pcs to a sparse matrix
    m = m.todok()  # Convert to Dictionary of Keys format for manipulation

    # Create Seurat-like object
    obj = {"counts": m, "meta_data": meta, "assay": "RNA", "graphs": {}}

    harmony_embeddings_all = pcs
    harmony_embeddings_all.index = range(harmony_embeddings_all.shape[0])

    # Create DimReducObject for harmony
    obj["reductions"] = {
        "harmony": {
            "embeddings": harmony_embeddings_all,
            "stdev": np.std(harmony_embeddings_all, axis=0, ddof=1),
            "assay": "RNA",
            "key": f"HARMONY",
        }
    }

    # Find neighbors (similar to Seurat)
    obj = find_neighbors(
        obj, reduction="harmony", dims=range(1, 21), k_param=30, nn_eps=0
    )

    # (1) Format data
    covs_keep = []
    if batches is not None:
        covs_keep.append(batches)
    if covs is not None:
        covs_keep.append(covs)

    if not obj["graphs"]:
        raise ValueError("Must precompute graph in Seurat with FindNeighbors()")

    if graph_use is None:
        graph_use = list(obj["graphs"].keys())[0]
        if verbose:
            print(f"Graph not specified. Using graph {graph_use}")
    else:
        if graph_use not in obj["graphs"]:
            raise ValueError(f"Graph {graph_use} not in seurat object")

    covs_keep.append(samplem_key)

    samplem_df = pd.DataFrame(
        obj["meta_data"][covs_keep].drop_duplicates()
    ).reset_index(drop=True)
    obs_df = obj["meta_data"].reset_index().rename(columns={"index": "CellID"})

    if samplem_df.shape[0] == obs_df.shape[0]:
        raise ValueError(
            "Sample-level metadata is the same length as cell-level metadata. "
            "Please check that samplem_vars are sample-level covariates."
        )

    rcna_data = {
        "samplem": samplem_df,
        "obs": obs_df,
        "connectivities": obj["graphs"][graph_use],
        "samplem_key": samplem_key,
        "obs_key": "CellID",
        "N": samplem_df.shape[0],
    }

    data = rcna_data
    suffix = ""

    # Formatting and error checking
    if batches is None:
        batches_vec = np.ones(data["N"], dtype=int)
    else:
        batches_vec = data["samplem"][batches].astype(int).values.flatten()

    max_frac_pcs = 0.15
    res = {}
    covs_mat = data["samplem"][covs].values

    # Formula generation (could be adapted based on specific use case)
    f = f"~0+{data['samplem_key']}"
    s = pd.get_dummies(
        data["obs"][data["obs_key"]].str.extract(f"^{data['samplem_key']}(.*)")[0]
    ).values
    s = s[data["samplem"][data["samplem_key"]].values]  # Necessary?

    # Define diffusion step function
    def diffuse_step(data, s):
        a = data["connectivities"]
        degrees = a.sum(axis=0) + 1
        s_norm = s / degrees  # Normalize
        res = (a @ s_norm) + s_norm  # Diffusion step
        return res

    prevmedkurt = np.inf
    maxnsteps = 15

    # Main diffusion loop
    for i in range(1, maxnsteps + 1):
        print(i)
        s = diffuse_step(data, s)
        print(s)

        medkurt = np.median(kurtosis(normalize(s, axis=1), axis=None))

        if nsteps is None:
            if prevmedkurt - medkurt < 3 and i > 3:
                print(f"stopping after {i} steps")
                break
            prevmedkurt = medkurt
            print(prevmedkurt)
            print(medkurt)
        elif i == nsteps:
            break

    return s  # or the relevant output


def find_neighbors(obj, reduction, dims, k_param, nn_eps):
    # Placeholder for finding neighbors
    # You will need to implement your own method to find neighbors based on the graph data
    obj["graphs"][f"{reduction}_neighbors"] = np.random.rand(
        len(obj["counts"]), len(obj["counts"])
    )  # Dummy placeholder
    return obj
