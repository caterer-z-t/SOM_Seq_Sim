import numpy as np


def tail_counts(z, znull):
    counts = []
    for znulli in znull.T:
        # Calculate counts for each znulli
        tail_counts = len(znulli) - np.cumsum(
            np.histogram(np.square(znulli), bins=np.append(0, np.square(z)))[0]
        )
        counts.append(tail_counts)
    return np.array(counts)
