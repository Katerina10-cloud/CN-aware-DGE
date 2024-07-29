from abc import ABC
from abc import abstractmethod
from typing import Literal
from typing import Optional
from typing import Tuple

import numpy as np
import pandas as pd


class Inference(ABC):
    """Abstract class with DESeq2-related inference methods."""
    
    @abstractmethod
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
        r"""Fit a NB GLM wit log-link to predict counts from the design matrix.

        See equations (1-2) in the DESeq2 paper.

        Parameters
        ----------
        counts : ndarray
            Raw counts.

        size_factors : ndarray
            Sample-wise scaling factors (obtained from median-of-ratios).

        design_matrix : ndarray
            Design matrix.

        disp : ndarray
            Gene-wise dispersion prior.

        min_mu : ndarray
            Lower bound on estimated means, to ensure numerical stability.
            (default: ``0.5``).

        beta_tol : float
            Stopping criterion for IRWLS:
            :math:`\vert dev - dev_{old}\vert / \vert dev + 0.1 \vert < \beta_{tol}`.
            (default: ``1e-8``).

        min_beta : float
            Lower-bound on LFC. (default: ``-30``).

        max_beta : float
            Upper-bound on LFC. (default: ``-30``).

        optimizer : str
            Optimizing method to use in case IRLS starts diverging.
            Accepted values: 'BFGS' or 'L-BFGS-B'.
            NB: only 'L-BFGS-B' ensures that LFCS will
            lay in the [min_beta, max_beta] range. (default: ``'L-BFGS-B'``).

        maxiter : int
            Maximum number of IRLS iterations to perform before switching to L-BFGS-B.
            (default: ``250``).

        Returns
        -------
        beta: ndarray
            Fitted (basemean, lfc) coefficients of negative binomial GLM.

        mu: ndarray
            Means estimated from size factors and beta:
            :math:`\mu = s_{ij} \exp(\beta^t X)`.

        H: ndarray
            Diagonal of the :math:`W^{1/2} X (X^t W X)^-1 X^t W^{1/2}`
            covariance matrix.

        converged: ndarray
            Whether IRLS or the optimizer converged. If not and if dimension allows it,
            perform grid search.
        """    

    @abstractmethod
    def alpha_mle(
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
        """Estimate the dispersion parameter of a negative binomial GLM.

        Parameters
        ----------
        counts : ndarray
            Raw counts.

        design_matrix : ndarray
            Design matrix.

        mu : ndarray
            Mean estimation for the NB model.

        alpha_hat : ndarray
            Initial dispersion estimate.

        min_disp : float
            Lower threshold for dispersion parameters.

        max_disp : float
            Upper threshold for dispersion parameters.

        prior_disp_var : float
            Prior dispersion variance.

        cr_reg : bool
            Whether to use Cox-Reid regularization. (default: ``True``).

        prior_reg : bool
            Whether to use prior log-residual regularization. (default: ``False``).

        optimizer : str
            Optimizing method to use. Accepted values: 'BFGS' or 'L-BFGS-B'.
            (default: ``'L-BFGS-B'``).

        Returns
        -------
        ndarray
            Dispersion estimate.

        ndarray
            Whether L-BFGS-B converged. If not, dispersion is estimated
            using grid search.
        """

    @abstractmethod
    def fit_rough_dispersions(
        self, normed_counts: np.ndarray, design_matrix: np.ndarray
    ) -> np.ndarray:
        """'Rough dispersion' estimates from linear model, as per the R code.

        Used as initial estimates in :meth:`DeseqDataSet.fit_genewise_dispersions()
        <pydeseq2.dds.DeseqDataSet.fit_genewise_dispersions>`.

        Parameters
        ----------
        normed_counts : ndarray
            Array of deseq2-normalized read counts. Rows: samples, columns: genes.

        design_matrix : pandas.DataFrame
            A DataFrame with experiment design information (to split cohorts).
            Indexed by sample barcodes. Unexpanded, *with* intercept.

        Returns
        -------
        ndarray
            Estimated dispersion parameter for each gene.
        """
    
    
    @abstractmethod 
    def fit_moments_dispersions(
        self, normed_counts: np.ndarray, size_factors: np.ndarray
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


