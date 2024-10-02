import numpy as np
from tail_counts import tail_counts


def empirical_fdrs(z, znull, thresholds):
    n = len(thresholds) - 1
    tails = tail_counts(thresholds, znull)[:n, :]
    ranks = tail_counts(thresholds, z)[:n, :]

    # Compute FDPs
    fdp = tails / ranks
    fdr = np.nanmean(
        fdp, axis=1
    )  # Handle division by zero by averaging only valid entries

    return fdr
