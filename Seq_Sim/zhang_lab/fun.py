import numpy as np


def fun(x, n, frac):
    min_by_frac = int(len(x) * frac)

    if len(x) <= n:
        return x
    elif len(x) > n and min_by_frac <= n:
        return np.random.choice(x, n, replace=False)
    elif len(x) > n and min_by_frac > n:
        return np.random.choice(x, min_by_frac, replace=False)
