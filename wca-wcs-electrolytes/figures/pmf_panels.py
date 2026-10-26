"""Loads a potential of mean force from one of the data/simulation/opes_pmf archives."""
import numpy as np


def load(path, prefix=None):
    z = np.load(path)
    if prefix:
        return z["cen"], z[f"{prefix}_mean"], z[f"{prefix}_sd"], z[f"{prefix}_cov"]
    return z["cen"], z["mean"], z["sd"], z["cov"]
