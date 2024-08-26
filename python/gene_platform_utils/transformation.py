"""
Variance stabilizing transformation (VST) for RNA-seq data like DESeq2's vst function
with blind=TRUE.

Adapted from the pydeseq2 package: https://github.com/owkin/PyDESeq2

MIT License

A pseudo-count value of 1 is added to the counts to avoid undefined log(0),
simplify operations and avoid resorting to a more costly iterative method.
Compare the output to DESeq2 on your data to check the difference.

"""

from typing import Optional, Union
import warnings

import numpy as np
from numpy.typing import NDArray
from scipy.optimize import minimize
from scipy.special import polygamma
from scipy.stats import trim_mean, norm
from scipy.special import gammaln


def check_dims(counts: np.ndarray) -> None:
    if counts.shape[0] > counts.shape[1]:
        warnings.warn(
            f"Found {counts.shape[1]} genes and {counts.shape[0]} samples."
            " If it's actually the other way around, you need to transpose the counts matrix"
            " before calling this function!",
            UserWarning,
        )


def vst(
    counts: np.ndarray,
    min_mu: float = 0.5,
    min_disp: float = 1e-8,
    max_disp: float = 10.0,
    beta_tol: float = 1e-8,
    subset_n_genes: Optional[int] = 1200,
    rand_seed: Optional[int] = 123,
) -> NDArray:
    """Variance stabilizing transformation (VST).

    Args:
        counts: The counts matrix.
        min_mu: The minimum mu.
        min_disp: The minimum dispersion.
        max_disp: The maximum dispersion.
        beta_tol: The beta tolerance.
        subset_n_genes: The number of genes to subset for faster dispersion fitting.
                        If None, all genes are used. The default is 1200.
        rand_seed: The random seed for gene subsetting. The default is 123.

    Returns:
        VST-transformed counts.

    """
    check_dims(counts)
    X = counts + 1
    max_disp = max(max_disp, len(X))

    normed_counts, size_factors, normed_means = fit_size_factors(X)

    if subset_n_genes is None:
        subset_n_genes = X.shape[1]
    else:
        subset_n_genes = min(X.shape[1], subset_n_genes)
    rng = np.random.default_rng(rand_seed)
    gene_subset = np.sort(rng.choice(X.shape[1], subset_n_genes, replace=False))

    MoM_dispersions = fit_MoM_dispersions(
        normed_counts[:, gene_subset], size_factors, min_disp, max_disp
    )

    dispersions = fit_dispersions(
        X[:, gene_subset],
        size_factors,
        min_mu,
        min_disp,
        max_disp,
        beta_tol,
        MoM_dispersions,
    )

    a0, a1 = fit_parametric_dispersion_trend(
        dispersions, normed_means[gene_subset], min_disp
    )

    return np.log2(
        (
            1
            + a1
            + 2 * a0 * normed_counts
            + 2 * np.sqrt(a0 * normed_counts * (1 + a1 + a0 * normed_counts))
        )
        / (4 * a0)
    )


def log2_1p(counts: np.ndarray) -> np.ndarray:
    return np.log2(counts + 1)


def ln_1p(counts: np.ndarray) -> np.ndarray:
    return np.log1p(counts)


def log10_1p(counts: np.ndarray) -> np.ndarray:
    return np.log10(counts + 1)


def fit_size_factors(X: NDArray) -> tuple[NDArray, NDArray, NDArray]:
    """Fit size factors.

    Args:
        X: The counts matrix.

    Returns:
        The normed counts, size factors, and normed means.

    """
    log_X = np.log(X)
    logmeans = log_X.mean(0)
    log_ratios = log_X - logmeans
    log_medians = np.median(log_ratios, axis=1)
    size_factors = np.exp(log_medians)
    normed_counts = X / size_factors[:, None]
    normed_means = normed_counts.mean(0)
    return normed_counts, size_factors, normed_means


def fit_MoM_dispersions(
    normed_counts: NDArray,
    size_factors: NDArray,
    min_disp: float,
    max_disp: float,
) -> NDArray:
    """Fit method of moments (MoM) dispersions.

    Args:
        normed_counts: The normalized counts.
        size_factors: The size factors.
        min_disp: The minimum dispersion.
        max_disp: The maximum dispersion.

    Returns:
        The MoM dispersions.

    """
    num_samples = normed_counts.shape[0]
    means = np.maximum(normed_counts.mean(0), 1)[None]
    alpha_rde = (
        ((normed_counts - means) ** 2 - means) / ((num_samples - 1) * means**2)
    ).sum(0)
    rde = np.maximum(alpha_rde, 0)
    s_mean_inv = (1 / size_factors).mean()
    mu = normed_counts.mean(0)
    sigma = normed_counts.var(0, ddof=1)
    mde = (sigma - s_mean_inv * mu) / mu**2
    alpha_hat = np.minimum(rde, mde)
    return np.clip(alpha_hat, min_disp, max_disp)


def fit_dispersions(
    X: NDArray,
    size_factors: NDArray,
    min_mu: float,
    min_disp: float,
    max_disp: float,
    beta_tol: float,
    MoM_dispersions: NDArray,
) -> NDArray:
    """Fit genewise dispersions.

    Args:
        X: The counts matrix.
        size_factors: The size factors.
        min_mu: The minimum mu.
        min_disp: The minimum dispersion.
        max_disp: The maximum dispersion.
        beta_tol: The beta tolerance.
        MoM_dispersions: The method of moments dispersions.

    Returns:
        The genewise dispersions.

    """
    _, mu_hat_, _, _ = irls(
        counts=X,
        size_factors=size_factors,
        disp=MoM_dispersions,
        min_mu=min_mu,
        beta_tol=beta_tol,
    )

    dispersions = alpha_mle(
        counts=X,
        mu=mu_hat_,
        alpha_hat=MoM_dispersions,
        min_disp=min_disp,
        max_disp=max_disp,
    )

    return np.clip(dispersions, min_disp, max_disp)


def fit_mean_dispersion_trend(
    dispersions: NDArray, min_disp: float
) -> tuple[float, NDArray]:
    """Fit a mean dispersion trend.

    Args:
        dispersions: The dispersions.
        min_disp: The minimum dispersion.

    Returns:
        The mean dispersion and the fitted dispersions.

    """
    mean_disp = trim_mean(
        dispersions[dispersions > 10 * min_disp], proportiontocut=0.001
    )
    fitted_dispersions = np.full(dispersions.shape[0], mean_disp)
    return mean_disp, fitted_dispersions


def fit_parametric_dispersion_trend(
    dispersions: NDArray,
    normed_means: NDArray,
    min_disp: float,
) -> NDArray:
    """Fit a parametric dispersion trend.

    Args:
        dispersions: The dispersions.
        normed_means: The normalized means.
        min_disp: The minimum dispersion.

    Returns:
        The trend coefficients and the fitted dispersions.

    """
    targets = dispersions
    covariates = 1 / normed_means
    mask = np.ones(dispersions.shape[0], dtype=bool)

    for i in range(len(targets)):
        if np.isinf(covariates[i]).any() or np.isnan(covariates[i]).any():
            mask[i] = False

    old_coeffs = np.array([0.1, 0.1])
    coeffs = np.array([1.0, 1.0])
    while (coeffs > 1e-10).all() and (
        np.log(np.abs(coeffs / old_coeffs)) ** 2
    ).sum() >= 1e-6:
        old_coeffs = coeffs
        coeffs, predictions, converged = dispersion_trend_gamma_glm(
            covariates[mask], targets[mask]
        )

        if not converged or (coeffs <= 1e-10).any():
            warnings.warn("Switching to a mean-based dispersion trend.", UserWarning)
            return fit_mean_dispersion_trend(dispersions, min_disp)

        # Filter out genes that are too far away from the curve before refitting
        pred_ratios = dispersions[mask] / predictions
        mask[np.where(mask)[0][(pred_ratios < 1e-4) | (pred_ratios >= 15)]] = False

    return coeffs


def grid_fit_alpha(
    counts: np.ndarray,
    mu: np.ndarray,
    alpha_hat: float,
    min_disp: float,
    max_disp: float,
    prior_disp_var: Optional[float] = None,
    cr_reg: bool = True,
    prior_reg: bool = False,
    grid_length: int = 100,
) -> float:
    min_log_alpha = np.log(min_disp)
    max_log_alpha = np.log(max_disp)
    grid = np.linspace(min_log_alpha, max_log_alpha, grid_length)

    def loss(log_alpha: np.ndarray) -> np.ndarray:
        alpha = np.exp(log_alpha)
        W = mu[:, None] / (1 + mu[:, None] * alpha)
        reg = 0
        if cr_reg:
            reg += 0.5 * np.log(W.sum(0))
        if prior_reg:
            reg += (np.log(alpha) - np.log(alpha_hat)) ** 2 / (2 * prior_disp_var)
        return vec_nb_nll(counts, mu, alpha) + reg

    ll_grid = loss(grid)

    min_idx = np.argmin(ll_grid)
    delta = grid[1] - grid[0]
    fine_grid = np.linspace(grid[min_idx] - delta, grid[min_idx] + delta, grid_length)

    ll_grid = loss(fine_grid)

    min_idx = np.argmin(ll_grid)
    log_alpha = fine_grid[min_idx]
    return log_alpha


def dnb_nll(counts: np.ndarray, mu: np.ndarray, alpha: float) -> float:
    alpha_neg1 = 1 / alpha
    ll_part = (
        alpha_neg1**2
        * (
            polygamma(0, alpha_neg1)
            - polygamma(0, counts + alpha_neg1)
            + np.log(1 + mu * alpha)
            + (counts - mu) / (mu + alpha_neg1)
        ).sum()
    )

    return -ll_part


def fit_alpha_mle(
    counts: np.ndarray,
    mu: np.ndarray,
    alpha_hat: float,
    min_disp: float,
    max_disp: float,
    prior_disp_var: Optional[float] = None,
    cr_reg: bool = True,
    prior_reg: bool = False,
) -> tuple[float, bool]:
    log_alpha_hat = np.log(alpha_hat)

    def loss(log_alpha: float) -> float:
        alpha = np.exp(log_alpha)
        reg = 0
        if cr_reg:
            W = mu / (1 + mu * alpha)
            reg += 0.5 * np.log(W.sum())
        if prior_reg:
            reg += (log_alpha - log_alpha_hat) ** 2 / (2 * prior_disp_var)
        return nb_nll(counts, mu, alpha) + reg

    def dloss(log_alpha: float) -> float:
        alpha = np.exp(log_alpha)
        reg_grad = 0
        if cr_reg:
            W = mu / (1 + mu * alpha)
            reg_grad += 0.5 * alpha * -(W**2).sum() / W.sum()

        if prior_reg:
            reg_grad += (log_alpha - log_alpha_hat) / prior_disp_var

        return alpha * dnb_nll(counts, mu, alpha) + reg_grad

    res = minimize(
        lambda x: loss(x[0]),
        x0=np.log(alpha_hat),
        jac=lambda x: dloss(x[0]),
        method="L-BFGS-B",
        bounds=[(np.log(min_disp), np.log(max_disp))],
    )

    if res.success:
        return np.exp(res.x[0])

    return np.exp(grid_fit_alpha(counts, mu, alpha_hat, min_disp, max_disp))


def vec_nb_nll(counts: np.ndarray, mu: np.ndarray, alpha: np.ndarray) -> np.ndarray:
    n = counts.shape[0]
    alpha_neg1 = 1 / alpha
    logbinom = (
        gammaln(counts[:, None] + alpha_neg1)
        - gammaln(counts + 1)[:, None]
        - gammaln(alpha_neg1)
    )

    if len(mu.shape) == 1:
        return n * alpha_neg1 * np.log(alpha) + (
            -logbinom
            + (counts[:, None] + alpha_neg1) * np.log(mu[:, None] + alpha_neg1)
            - (counts * np.log(mu))[:, None]
        ).sum(0)
    else:
        return n * alpha_neg1 * np.log(alpha) + (
            -logbinom
            + (counts[:, None] + alpha_neg1) * np.log(mu + alpha_neg1)
            - (counts[:, None] * np.log(mu))
        ).sum(0)


def irls_solver(
    counts: np.ndarray,
    size_factors: np.ndarray,
    disp: float,
    min_mu: float = 0.5,
    beta_tol: float = 1e-8,
    min_beta: float = -30,
    max_beta: float = 30,
    maxiter: int = 250,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, bool]:
    y = np.log(counts / size_factors + 0.1)
    beta_init = y.sum()[None] / len(counts)
    beta = beta_init

    dev = 1000.0
    dev_ratio = 1.0

    ridge_factor = np.array([[1e-6]])
    mu = np.maximum(size_factors * np.exp(beta), min_mu)

    converged = True
    i = 0
    while dev_ratio > beta_tol:
        W = mu / (1.0 + mu * disp)
        z = np.log(mu / size_factors) + (counts - mu) / mu
        H = W[None].sum()[None] + ridge_factor
        beta_hat = ((W * z).sum()[None] / H)[0]
        i += 1

        if sum(np.abs(beta_hat) > max_beta) > 0 or i >= maxiter:

            def f(beta: np.ndarray) -> float:
                # closure to minimize
                mu_ = np.maximum(size_factors * np.exp(beta), min_mu)
                return nb_nll(counts, mu_, disp) + 0.5 * (ridge_factor @ beta**2).sum()

            def df(beta: np.ndarray) -> np.ndarray:
                mu_ = np.maximum(size_factors * np.exp(beta), min_mu)
                return (
                    -counts.sum()[None]
                    + ((1 / disp + counts) * mu_ / (1 / disp + mu_)).sum()[None]
                    + ridge_factor @ beta
                )

            res = minimize(
                f,
                beta_init,
                jac=df,
                method="L-BFGS-B",
                bounds=[(min_beta, max_beta)],
            )

            beta = res.x
            mu = np.maximum(size_factors * beta, min_mu)
            converged = res.success
            break

        beta = beta_hat
        mu = np.maximum(size_factors * np.exp(beta), min_mu)
        old_dev = dev
        dev = -2 * nb_nll(counts, mu, disp)
        dev_ratio = np.abs(dev - old_dev) / (np.abs(dev) + 0.1)

    W = mu / (1.0 + mu * disp)
    H = np.full_like(W, 1 / (np.sum(W) + 1e-6))
    W_sq = np.sqrt(W)
    H = W_sq * H * W_sq

    mu = size_factors * np.exp(beta)
    return beta, mu, H, converged


def irls(
    counts: np.ndarray,
    size_factors: np.ndarray,
    disp: np.ndarray,
    min_mu: float,
    beta_tol: float,
    min_beta: float = -30,
    max_beta: float = 30,
    maxiter: int = 250,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    res = (
        irls_solver(
            counts=counts[:, i],
            size_factors=size_factors,
            disp=disp[i],
            min_mu=min_mu,
            beta_tol=beta_tol,
            min_beta=min_beta,
            max_beta=max_beta,
            maxiter=maxiter,
        )
        for i in range(counts.shape[1])
    )

    res = zip(*res)
    MLE_lfcs_, mu_hat_, hat_diagonals_, converged_ = (np.array(m) for m in res)
    return (
        MLE_lfcs_,
        mu_hat_.T,
        hat_diagonals_.T,
        converged_,
    )


def alpha_mle(
    counts: np.ndarray,
    mu: np.ndarray,
    alpha_hat: np.ndarray,
    min_disp: float,
    max_disp: float,
    prior_disp_var: Optional[float] = None,
    cr_reg: bool = True,
    prior_reg: bool = False,
) -> tuple[np.ndarray, np.ndarray]:
    return np.array(
        [
            fit_alpha_mle(
                counts=counts[:, i],
                mu=mu[:, i],
                alpha_hat=alpha_hat[i],
                min_disp=min_disp,
                max_disp=max_disp,
                prior_disp_var=prior_disp_var,
                cr_reg=cr_reg,
                prior_reg=prior_reg,
            )
            for i in range(counts.shape[1])
        ]
    )


def mean_absolute_deviation(x: np.ndarray) -> float:
    center = np.median(x)
    return np.median(np.abs(x - center)) / norm.ppf(0.75)


def nb_nll(
    counts: np.ndarray, mu: np.ndarray, alpha: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    alpha_neg1 = 1 / alpha
    logbinom = gammaln(counts + alpha_neg1) - gammaln(counts + 1) - gammaln(alpha_neg1)
    return (
        counts.shape[0] * alpha_neg1 * np.log(alpha)
        + (
            -logbinom
            + (counts + alpha_neg1) * np.log(alpha_neg1 + mu)
            - counts * np.log(mu)
        ).sum()
    )


def dispersion_trend_gamma_glm(
    covariates: NDArray, targets: NDArray
) -> tuple[np.ndarray, np.ndarray, bool]:
    covariates = np.column_stack((np.ones(covariates.shape[0]), covariates))

    def loss(coeffs):
        mu = covariates @ coeffs
        return np.nanmean(targets / mu + np.log(mu), axis=0)

    def grad(coeffs):
        mu = covariates @ coeffs
        return -np.nanmean(
            ((targets / mu - 1)[:, None] * covariates) / mu[:, None], axis=0
        )

    try:
        res = minimize(
            loss,
            x0=np.array([1.0, 1.0]),
            jac=grad,
            method="L-BFGS-B",
            bounds=[(1e-12, np.inf)],
        )
    except RuntimeWarning:
        return np.array([np.nan, np.nan]), np.array([np.nan, np.nan]), False

    coeffs = res.x
    return coeffs, covariates @ coeffs, res.success
