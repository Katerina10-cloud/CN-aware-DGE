from typing import Literal
from typing import Optional
from typing import Tuple

import numpy as np
import pandas as pd
from joblib import Parallel  # type: ignore
from joblib import delayed
from joblib import parallel_backend
from scipy.optimize import minimize  # type: ignore

from pydeseq2 import utils
import inference
import utils_pydeseq2CN 

class DefInference(inference.Inference):

    """Default DESeq2-related inference methods, using scipy/sklearn/numpy.

    This object contains the interface to the default inference routines and uses
    joblib internally for parallelization. Inherit this class or its parent to write
    custom inference routines.

    Parameters
    ----------
    joblib_verbosity : int
        The verbosity level for joblib tasks. The higher the value, the more updates
        are reported. (default: ``0``).
    batch_size : int
        Number of tasks to allocate to each joblib parallel worker. (default: ``128``).
    n_cpus : int
        Number of cpus to use. If None, all available cpus will be used.
        (default: ``None``).
    backend : str
        Joblib backend.
    """
    
    fit_rough_dispersions = staticmethod(utils.fit_rough_dispersions)  # type: ignore
    fit_moments_dispersions = staticmethod(utils.fit_moments_dispersions)  # type: ignore
    
    def __init__(
        self,
        joblib_verbosity: int = 0,
        batch_size: int = 128,
        n_cpus: Optional[int] = None,
        backend: str = "loky",
    ):
        self._joblib_verbosity = joblib_verbosity
        self._batch_size = batch_size
        self._n_cpus = utils.get_num_processes(n_cpus)
        self._backend = backend
        
    @property
    def n_cpus(self) -> int:  # noqa: D102
        return self._n_cpus

    @n_cpus.setter
    def n_cpus(self, n_cpus: int) -> None:
        self._n_cpus = utils.get_num_processes(n_cpus)

    def irls_glm( 
        self,
        counts: np.ndarray,
        size_factors: np.ndarray,
        design_matrix: np.ndarray,
        disp: np.ndarray,
        cnv: np.ndarray,
        min_mu: float,
        beta_tol: float,
        min_beta: float = -30,
        max_beta: float = 30,
        optimizer: Literal["BFGS", "L-BFGS-B"] = "L-BFGS-B",
        maxiter: int = 250,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        with parallel_backend(self._backend, inner_max_num_threads=1):
            res = Parallel(
                n_jobs=self.n_cpus,
                verbose=self._joblib_verbosity,
                batch_size=self._batch_size,
            )(
                delayed(utils_pydeseq2CN.irls_glm)(
                    counts=counts[:, i],
                    size_factors=size_factors,
                    design_matrix=design_matrix,
                    disp=disp[i],
                    cnv=cnv[:, i],
                    min_mu=min_mu,
                    beta_tol=beta_tol,
                    min_beta=min_beta,
                    max_beta=max_beta,
                    optimizer=optimizer,
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
    def alpha_mle(  # noqa: D102
        self,
        counts: np.ndarray,
        design_matrix: np.ndarray,
        mu: np.ndarray,
        alpha_hat: np.ndarray,
        min_disp: float,
        max_disp: float,
        prior_disp_var: Optional[float] = None,
        cr_reg: bool = True,
        prior_reg: bool = False,
        optimizer: Literal["BFGS", "L-BFGS-B"] = "L-BFGS-B",
    ) -> Tuple[np.ndarray, np.ndarray]:
        with parallel_backend(self._backend, inner_max_num_threads=1):
            res = Parallel(
                n_jobs=self.n_cpus,
                verbose=self._joblib_verbosity,
                batch_size=self._batch_size,
            )(
                delayed(utils.fit_alpha_mle)(
                    counts=counts[:, i],
                    design_matrix=design_matrix,
                    mu=mu[:, i],
                    alpha_hat=alpha_hat[i],
                    min_disp=min_disp,
                    max_disp=max_disp,
                    prior_disp_var=prior_disp_var,
                    cr_reg=cr_reg,
                    prior_reg=prior_reg,
                    optimizer=optimizer,
                )
                for i in range(counts.shape[1])
            )
        res = zip(*res)
        dispersions_, l_bfgs_b_converged_ = (np.array(m) for m in res)
        return dispersions_, l_bfgs_b_converged_
