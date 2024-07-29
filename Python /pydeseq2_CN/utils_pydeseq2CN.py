import multiprocessing
import warnings
from math import ceil
from math import floor
from pathlib import Path
from typing import List
from typing import Literal
from typing import Optional
from typing import Tuple
from typing import Union
from typing import cast

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.linalg import solve  # type: ignore
from scipy.optimize import minimize  # type: ignore
from scipy.special import gammaln  # type: ignore
from scipy.special import polygamma  # type: ignore
from scipy.stats import norm  # type: ignore
from sklearn.linear_model import LinearRegression  # type: ignore

import pydeseq2
import pydeseq2.utils
from pydeseq2.grid_search import grid_fit_alpha
from pydeseq2.grid_search import grid_fit_shrink_beta

def irls_glm(
    counts: np.ndarray,
    size_factors: np.ndarray,
    design_matrix: np.ndarray,
    cnv: np.ndarray,
    disp: float,
    min_mu: float = 0.5,
    beta_tol: float = 1e-8,
    min_beta: float = -30,
    max_beta: float = 30,
    optimizer: Literal["BFGS", "L-BFGS-B"] = "L-BFGS-B",
    maxiter: int = 250,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, bool]:

    assert optimizer in ["BFGS", "L-BFGS-B"]
    
    X = design_matrix
    num_vars = design_matrix.shape[1]
    
    # if full rank, estimate initial betas for IRLS below
    if np.linalg.matrix_rank(X) == num_vars:
        Q, R = np.linalg.qr(X)
        y = np.log((counts/cnv)/size_factors + 0.1)
        beta_init = solve(R, Q.T @ y)
        beta = beta_init

    else:  # Initialise intercept with log base mean
        beta_init = np.zeros(num_vars)
        beta_init[0] = np.log((counts / cnv) / size_factors).mean()
        beta = beta_init
        
    dev = 1000.0
    dev_ratio = 1.0

    ridge_factor = np.diag(np.repeat(1e-6, num_vars))
    mu = np.maximum(cnv * size_factors * np.exp(X @ beta), min_mu)
    
    converged = True
    i = 0
    while dev_ratio > beta_tol:
        W = mu / (1.0 + mu * disp)
        z = np.log((mu / cnv)/size_factors) + (counts - mu) / mu
        H = (X.T * W) @ X + ridge_factor
        beta_hat = solve(H, X.T @ (W * z), assume_a="pos")
        i += 1

        if sum(np.abs(beta_hat) > max_beta) > 0 or i >= maxiter:
            # If IRLS starts diverging, use L-BFGS-B
            def f(beta: np.ndarray) -> float:
                # closure to minimize
                mu_ = np.maximum(cnv * size_factors * np.exp(X @ beta), min_mu)
                
                return nb_nll(counts, mu_, disp) + 0.5 * (ridge_factor @ beta**2).sum()

            def df(beta: np.ndarray) -> np.ndarray:
                mu_ = np.maximum(cnv * size_factors * np.exp(X @ beta), min_mu)
                #mu_ = np.maximum(size_factors * np.exp(X @ beta), min_mu)
                return (
                    -X.T @ counts
                    + ((1 / disp + counts) * mu_ / (1 / disp + mu_)) @ X
                    + ridge_factor @ beta
                )

            res = minimize(
                f,
                beta_init,
                jac=df,
                method=optimizer,
                bounds=(
                    [(min_beta, max_beta)] * num_vars
                    if optimizer == "L-BFGS-B"
                    else None
                ),
            )
            beta = res.x
            mu = np.maximum(cnv * size_factors * np.exp(X @ beta), min_mu)
            converged = res.success

            if not res.success and num_vars <= 2:
                beta = grid_fit_beta(
                    counts,
                    size_factors,
                    cnv,
                    X,
                    disp,
                )
                mu = np.maximum(cnv * size_factors * np.exp(X @ beta), min_mu)
            break

        beta = beta_hat
        mu = np.maximum(cnv * size_factors * np.exp(X @ beta), min_mu)
        
        # Compute deviation
        old_dev = dev
        # Replaced deviation with -2 * nll, as in the R code
        dev = -2 * nb_nll(counts, mu, disp)
        dev_ratio = np.abs(dev - old_dev) / (np.abs(dev) + 0.1)

    # Compute H diagonal (useful for Cook distance outlier filtering)
    W = mu / (1.0 + mu * disp)
    W_sq = np.sqrt(W)
    XtWX = (X.T * W) @ X + ridge_factor
    H = W_sq * np.diag(X @ np.linalg.inv(XtWX) @ X.T) * W_sq
    
    # Return an UNthresholded mu (as in the R code)
    # Previous quantities are estimated with a threshold though
    
    return beta, mu, H, converged



def fit_moments_dispersions(
    normed_counts: np.ndarray, size_factors: np.ndarray, cnv: np.ndarray
) -> np.ndarray:
    """Dispersion estimates based on moments, as per the R code.

    Used as initial estimates in :meth:`DeseqDataSet.fit_genewise_dispersions()
    <pydeseq2.dds.DeseqDataSet.fit_genewise_dispersions>`.

    Parameters
    ----------
    normed_counts : ndarray
        Array of deseq2-normalized read counts. Rows: samples, columns: genes.

    size_factors : ndarray
        DESeq2 normalization factors.

    Returns
    -------
    ndarray
        Estimated dispersion parameter for each gene.
    """
    # Exclude genes with all zeroes
    normed_counts = normed_counts[:, ~(normed_counts == 0).all(axis=0)]
    # mean inverse size factor
    s_mean_inv = ((1 / cnv)/size_factors).mean()
    mu = normed_counts.mean(0)
    sigma = normed_counts.var(0, ddof=1)
    # ddof=1 is to use an unbiased estimator, as in R
    # NaN (variance = 0) are replaced with 0s
    return np.nan_to_num((sigma - s_mean_inv * mu) / mu**2)



def grid_fit_beta(
    counts: np.ndarray,
    size_factors: np.ndarray,
    design_matrix: np.ndarray,
    disp: float,
    cnv: np.ndarray,
    min_mu: float = 0.5,
    grid_length: int = 60,
    min_beta: float = -30,
    max_beta: float = 30,
) -> np.ndarray:
    
    x_grid = np.linspace(min_beta, max_beta, grid_length)
    y_grid = np.linspace(min_beta, max_beta, grid_length)
    ll_grid = np.zeros((grid_length, grid_length))

    def loss(beta: np.ndarray) -> np.ndarray:
        # closure to minimize

        mu = np.maximum(cnv[:, None] * size_factors[:, None] * np.exp(design_matrix @ beta.T), min_mu)
        return vec_nb_nll(counts, mu, disp) + 0.5 * (1e-6 * beta**2).sum(1)

    for i, x in enumerate(x_grid):
        ll_grid[i, :] = loss(np.array([[x, y] for y in y_grid]))

    min_idxs = np.unravel_index(np.argmin(ll_grid, axis=None), ll_grid.shape)
    delta = x_grid[1] - x_grid[0]

    fine_x_grid = np.linspace(
        x_grid[min_idxs[0]] - delta, x_grid[min_idxs[0]] + delta, grid_length
    )

    fine_y_grid = np.linspace(
        y_grid[min_idxs[1]] - delta,
        y_grid[min_idxs[1]] + delta,
        grid_length,
    )

    for i, x in enumerate(fine_x_grid):
        ll_grid[i, :] = loss(np.array([[x, y] for y in fine_y_grid]))

    min_idxs = np.unravel_index(np.argmin(ll_grid, axis=None), ll_grid.shape)
    beta = np.array([fine_x_grid[min_idxs[0]], fine_y_grid[min_idxs[1]]])
    return beta
