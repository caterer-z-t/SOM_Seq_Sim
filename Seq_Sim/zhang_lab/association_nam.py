import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import kurtosis
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import TruncatedSVD
import random
import warnings


def association_nam(
    seurat_object=None,
    metadata=None,
    pcs=None,
    test_var=None,
    samplem_key=None,
    graph_use="RNA_snn",
    batches=None,
    covs=None,
    nsteps=None,
    verbose=True,
    assay=None,
    key="NAMPC_",
    maxnsteps=15,
    max_frac_pcs=0.15,
    ks=None,
    Nnull=1000,
    force_permute_all=False,
    local_test=True,
    seed=1234,
    return_nam=True,
):
    if seurat_object is not None and metadata is None and pcs is None:
        print("Will use Seurat object following analysis...")
    elif seurat_object is None and metadata is not None and pcs is not None:
        # Create a Seurat object equivalent in Python
        meta = metadata.copy()
        meta.index = range(len(meta))

        m = sparse.csr_matrix(pcs.T)  # Convert pcs to a sparse matrix
        obj = {
            "counts": m,
            "meta_data": meta,
            "assay": "RNA",
        }

        harmony_embeddings_all = pcs
        obj["reductions"] = {
            "harmony": {
                "embeddings": harmony_embeddings_all,
                "stdev": np.std(harmony_embeddings_all, axis=0),
                "assay": "RNA",
                "key": "HARMONY",
            }
        }

        # Simulate FindNeighbors
        seurat_object = find_neighbors(
            obj, reduction="harmony", dims=range(20), k_param=30
        )

    elif (seurat_object is None and metadata is None and pcs is not None) or (
        seurat_object is None and metadata is not None and pcs is None
    ):
        raise ValueError("Must provide both metadata and precomputed PCs")

    # (1) Format data
    covs_keep = [test_var]
    if batches is not None:
        covs_keep.append(batches)
    if covs is not None:
        covs_keep.extend(covs)

    if not seurat_object.get("graphs"):
        raise ValueError("Must precompute graph in Seurat with FindNeighbors()")

    if graph_use not in seurat_object["graphs"]:
        raise ValueError(f"Graph {graph_use} not in seurat object")

    covs_keep.append(samplem_key)
    samplem_df = (
        pd.DataFrame(seurat_object["meta_data"][covs_keep])
        .drop_duplicates()
        .reset_index(drop=True)
    )
    obs_df = seurat_object["meta_data"].reset_index()

    if len(samplem_df) == len(obs_df):
        raise ValueError(
            "Sample-level metadata is the same length as cell-level metadata. Please check that samplem_vars are sample-level covariates."
        )

    rcna_data = {
        "samplem": samplem_df,
        "obs": obs_df,
        "connectivities": seurat_object["graphs"][graph_use],
        "samplem_key": samplem_key,
        "obs_key": "CellID",
        "N": len(samplem_df),
    }

    # Formatting and error checking
    if batches is None:
        batches_vec = np.ones(data["N"], dtype=int)
    else:
        batches_vec = np.array(seurat_object["meta_data"][batches]).astype(int)

    res = {}
    covs_mat = np.array(seurat_object["meta_data"][covs])

    s = pd.get_dummies(
        obs_df[samplem_key], prefix="", prefix_sep=""
    ).values  # Model matrix
    s = s[:, samplem_df[samplem_key].values]

    prevmedkurt = float("inf")

    def diffuse_step(data, s):
        a = data["connectivities"]
        degrees = np.array(a.sum(axis=1)).flatten() + 1
        s_norm = s / degrees[:, np.newaxis]
        return (a @ s_norm) + s_norm

    for i in range(maxnsteps):
        s = diffuse_step(rcna_data, s)
        medkurt = np.median(kurtosis(s, axis=1))
        if nsteps is None:
            prevmedkurt = medkurt
            if prevmedkurt - medkurt < 3 and i > 3:
                if verbose:
                    print(f"Stopping after {i} steps")
                break
        elif i == nsteps:
            break

    snorm = s / s.sum(axis=1, keepdims=True)  # Normalization
    NAM = snorm

    N = NAM.shape[0]
    if len(np.unique(batches_vec)) == 1:
        if verbose:
            print("Only one unique batch supplied to qc")
        keep = np.arange(N)
    else:
        if verbose:
            print("Filtering based on batches kurtosis")

    def batch_kurtosis(NAM, batches_vec):
        return pd.Series(
            [
                kurtosis(NAM[batches_vec == b].mean(axis=0))
                for b in np.unique(batches_vec)
            ]
        )

    kurtoses = batch_kurtosis(NAM, batches_vec)
    threshold = max(6, 2 * np.median(kurtoses))
    if verbose:
        print(f"Throwing out neighborhoods with batch kurtosis >= {threshold}")
    keep = np.where(kurtoses < threshold)[0]

    # Prepare results
    res[f"NAM.T"] = NAM.T
    res[f"keptcells"] = keep
    res[f"_batches"] = batches_vec

    if verbose:
        print("Residualize NAM")

    NAM_ = StandardScaler(with_mean=True, with_std=False).fit_transform(NAM)
    ncols_C = 0

    if covs_mat is not None:
        covs_mat = StandardScaler().fit_transform(covs_mat)
        ncols_C += covs_mat.shape[1]

    if covs_mat is None:
        M = sparse.eye(N)
    else:
        M = (
            sparse.eye(N)
            - covs_mat @ np.linalg.pinv(covs_mat.T @ covs_mat) @ covs_mat.T
        )

    NAM_ = M @ NAM_
    res[f"_M"] = M
    res[f"_r"] = ncols_C

    if verbose:
        print("Decompose NAM")

    npcs = min(10, round(max_frac_pcs * len(rcna_data["samplem"])))

    if npcs > 0.5 * min(NAM_.shape):
        svd = TruncatedSVD(n_components=npcs)
        svd_res = svd.fit(NAM_)
    else:
        svd_res = TruncatedSVD(n_components=npcs).fit(NAM_)

    U_df = svd_res.transform(NAM_)
    V_df = svd_res.components_.T

    res[f"NAM_sampleXpc"] = U_df
    res[f"NAM_svs"] = svd_res.singular_values_**2
    res[f"NAM_varexp"] = res[f"NAM_svs"] / len(U_df) / len(V_df)
    res[f"NAM_nbhdXpc"] = V_df

    NAMsvd = {
        "sampleXpc": res[f"NAM_sampleXpc"],
        "svs": res[f"NAM_svs"],
        "nbhdXpc": res[f"NAM_nbhdXpc"],
        "varexp": res[f"NAM_varexp"],
    }

    M = res[f"_M"]
    r = res[f"_r"]

    yvals = rcna_data["samplem"][test_var].values
    if pd.api.types.is_numeric_dtype(yvals):
        y = yvals
    else:
        raise ValueError(
            f"test_var is of class {type(yvals)}. It must be numeric variable for association testing."
        )

    random.seed(seed)
    if force_permute_all:
        batches_vec = np.ones(len(y), dtype=int)

    # Prepare U and other necessary components...
    # More code would be here for final processing of U

    return NAMsvd  # Modify according to your needs
