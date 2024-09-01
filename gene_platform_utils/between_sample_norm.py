"""edgeR-equivalent TMM normalization method adapted from github.com/genialis/RNAnorm"""

import numpy as np
from numpy.typing import NDArray
from scipy.stats import gmean, rankdata, scoreatpercentile


def compute_tmm_effective_library_sizes(
    X: NDArray[np.integer], m_trim: float = 0.3, a_trim: float = 0.05
) -> NDArray[np.floating]:
    """Get TMM normalization factors (un-normalized with geometric mean).

    Args:
        X:  Expression raw count matrix (n_samples, n_features)
        m_trim: Keep genes within (m_trim, 1 - m_trim) percentile of M-values.
        a_trim: Keep genes within (a_trim, 1 - a_trim) percentile of A-values.

    """
    # Remove genes with zero count in all samples and convert to float
    X = X[:, np.sum(X, axis=0) > 0].astype(np.float64)

    lib_size = np.nansum(X, axis=1)

    # Get reference sample using upper quartile norm factors
    upper_quartiles = np.apply_along_axis(scoreatpercentile, axis=1, arr=X, per=75)
    uq_factors = upper_quartiles / lib_size
    f75 = uq_factors / gmean(uq_factors)
    ref = X[np.argmin(np.fabs(f75 - np.mean(f75))), :]

    # Compute the effective library size of the reference sample
    lib_size_ref = np.nansum(ref[np.newaxis, :], axis=1)

    X[X == 0] = np.nan
    ref[ref == 0] = np.nan

    r = X / lib_size[:, np.newaxis]
    r_ref = ref / lib_size_ref

    M = np.log2(r / r_ref)
    A = np.log2(r * r_ref) / 2
    W = (1 - r) / X + (1 - r_ref) / ref

    finite = np.isfinite(M) & np.isfinite(A)
    N = np.sum(finite, axis=1)
    M_low = np.floor(N * m_trim) + 1
    M_high = N - M_low + 1
    A_low = np.floor(N * a_trim) + 1
    A_high = N - A_low + 1

    f = []
    for i in range(X.shape[0]):  # for each sample
        m, a, w = [M[i, finite[i]], A[i, finite[i]], W[i, finite[i]]]
        m_rank, a_rank = rankdata(m), rankdata(a)
        keep_genes_mask = np.logical_and(
            np.logical_and(m_rank >= M_low[i], m_rank <= M_high[i]),
            np.logical_and(a_rank >= A_low[i], a_rank <= A_high[i]),
        )
        f.append(np.nansum(keep_genes_mask * m / w) / np.nansum(keep_genes_mask / w))

    tmm_norm_factors = np.power(2, f)

    return np.nansum(X, axis=1) * tmm_norm_factors / gmean(tmm_norm_factors)
