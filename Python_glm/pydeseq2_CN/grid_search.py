from typing import Optional

import numpy as np
from scipy.special import gammaln  # type: ignore

import utils_pydeseq2CN


def grid_fit_shrink_beta(
    counts: np.ndarray,
    cnv: np.ndarray,
    offset: np.ndarray,
    design_matrix: np.ndarray,
    size: np.ndarray,
    prior_no_shrink_scale: float,
    prior_scale: float,
    scale_cnst: float,
    grid_length: int = 60,
    min_beta: float = -30,
    max_beta: float = 30,
) -> np.ndarray:
    """Find best LFC parameter.

    Performs 2D grid search to maximize MAP negative binomial
    GLM log-likelihood w.r.t. LFCs, with apeGLM prior.

    Parameters
    ----------
    counts : ndarray
        Raw counts for a given gene.

    offset : ndarray
        Natural logarithm of size factor.

    design_matrix : ndarray
        Design matrix.

    size : ndarray
        Size parameter of NB family (inverse of dispersion).

    prior_no_shrink_scale : float
        Prior variance for the intercept.

    prior_scale : float
        Prior variance for the LFC coefficient.

    scale_cnst : float
        Scaling factor for the optimization.

    grid_length : int
        Number of grid points. (default: ``100``).

    min_beta : int
        Lower-bound on LFC. (default: ``30``).

    max_beta : int
        Upper-bound on LFC. (default: ``30``).

    Returns
    -------
    ndarray
        Fitted MAP LFC parameter.
    """
    x_grid = np.linspace(min_beta, max_beta, grid_length)
    y_grid = np.linspace(min_beta, max_beta, grid_length)
    ll_grid = np.zeros((grid_length, grid_length))

    def loss(beta: np.ndarray) -> float:
        # closure to minimize
        return (
            utils_pydeseq2CN.nbinomFn(
                beta,
                design_matrix,
                counts,
                cnv,
                size,
                offset,
                prior_no_shrink_scale,
                prior_scale,
            )
            / scale_cnst
        )

    for i, x in enumerate(x_grid):
        for j, y in enumerate(y_grid):
            ll_grid[i, j] = loss(np.array([x, y]))

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
        for j, y in enumerate(fine_y_grid):
            ll_grid[i, j] = loss(np.array([x, y]))

    min_idxs = np.unravel_index(np.argmin(ll_grid, axis=None), ll_grid.shape)
    beta = np.array([fine_x_grid[min_idxs[0]], fine_y_grid[min_idxs[1]]])
    return beta