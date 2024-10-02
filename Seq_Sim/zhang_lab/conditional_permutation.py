import numpy as np
import pandas as pd


def conditional_permutation(B, Y, num):
    result = []
    for i in range(num):
        # Split indices based on B
        split_indices = {b: np.where(B == b)[0] for b in np.unique(B)}
        permuted_vals = []

        for idx in split_indices.values():
            # Sample values from Y based on the indices
            permuted_vals.append(np.random.permutation(Y[idx]))

        # Combine the results into a DataFrame
        combined = pd.DataFrame(
            {
                "idx": np.concatenate(list(split_indices.values())),
                "val": np.concatenate(permuted_vals),
            }
        )
        combined = combined.sort_values(by="idx")
        result.append(combined["val"].to_numpy())

    return np.column_stack(result)
