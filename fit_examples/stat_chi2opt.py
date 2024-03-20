r"""
    Author: Tong Li, Mich State Univ, 2021-
"""

import numpy as np
from copy import deepcopy
import functools
printf = functools.partial(print, flush=True) # force output flush
try:
    from mpi4py import MPI
    import_MPI = True
except ImportError:
    import_MPI = False


def _unpack_tuple(tin, n_param):
    assert(len(tin)==2)
    lb, ub = [np.asarray(b, dtype=float) for b in tin]
    if lb.ndim == 0: lb = np.resize(lb, n_param)
    if ub.ndim == 0: ub = np.resize(ub, n_param)
    assert((lb.ndim == 1) and (ub.ndim == 1) and (len(lb) == len(ub) == n_param))
    return lb, ub


def stat_from_derivatives(residual_2norm, J, outfile=None):
    r"""Calculate covariance, correlation and sensitivity matrix at the minimum point
    of $$\chi^2$$ using J, derivatives of residuals w.r.t. fitting parameters.
    J has the dimension of (n_param, n_obs), where n_param is the number of fitting parameters
    and n_obs is the number of observables / model predictions / residuals. 

    Parameters
    ----------
    residual_2norm: float
        2-norm of the residual vector, i.e. the sum of squared residuals.
        It is used to compute $$\chi^2$$, an estimate of the variance.
    J: array-like with shape (n_param, n_obs)
        Derivatives of residuals w.r.t. fitting parameters.
    outfile: IO object or None, optional.
        IO object to write output. If None (default), no output will be provided.
    
    Returns
    -------
    cov: ndarray with shape (n_param, n_param)
        Covariance matrix.
    eigvals: ndarray with shape (n_param,)
        Eigenvalues of the covariance matrix. 
    eigvecs: ndarray with shape (n_param, n_param)
        Eigenvectors of the covariance matrix. 
    corr : ndarray with shape (n_param, n_param)
        Correlation matrix.
    sen: ndarray with shape (n_param, n_obs)
        Sensitivity matrix. 
    
    References
    ----------
    .. M. Kortelainen, T. Lesinski, J. Moré, W. Nazarewicz, J. Sarich, N. Schunck, M. V. Stoitsov, and S. Wild,
       "Nuclear energy density optimization", Phys. Rev. C, 82, 024313 (2010).
    .. J. R. Donaldson and R. B. Schnabel, "Computational Experience With Confidence Regions and Confidence 
       Intervals for Nonlinear Least Squares", Technometrics 29, 67 (1987).
    """

    J = np.asarray(J, dtype=float)
    assert(J.ndim == 2)
    n_param, n_obs = J.shape
    dof = n_obs - n_param # degrees of freedom
    chi_square_min = residual_2norm / dof # minimum chi square
    cov = chi_square_min * np.linalg.inv(np.matmul(J, J.T)) # covariance matrix
    eigvals, eigvecs = np.linalg.eigh(cov) # eigenvalue decomposition of covariance
    corr = np.zeros(cov.shape) # correlation matrix
    for i in range(corr.shape[0]):
        for j in range(corr.shape[1]):
            if i==j:
                corr[i][j] = 1.0
            else:
                corr[i][j] = cov[i][j] / np.sqrt(cov[i][i]*cov[j][j])
    sensitivity = np.matmul(np.linalg.inv(np.matmul(J, J.T)), J) # sensitivity matrix
    if outfile is not None:
        printf('# J matrix: \n', J.tolist(), file=outfile)
        printf('# Covariance matrix: \n', cov.tolist(), file=outfile)
        printf('# Eigenvalues of covariance matrix: \n', eigvals.tolist(), file=outfile)
        printf('# Eigenvectors of covariance matrix: \n', eigvecs.tolist(), file=outfile)
        printf('# Correlation matrix: \n', corr.tolist(), file=outfile)
        printf('# Sensitivity matrix: \n', sensitivity.tolist(), file=outfile)
    return cov, eigvals, eigvecs, corr, sensitivity


def stat_analysis(fun, x, diff, res_min=None, bounds=(-np.inf, np.inf), 
                  output_file='stat_analysis.txt', args=(), kwargs={}, comm=None):
    r"""Use the finite difference method to compute derivatives around the 
    minimum point obtained from POUNDerS or other least square minimization routines. 
    Then the covariance, correlation and sensitivity matrix are calculated, based on
    equations given in Phys. Rev. C 82 024313 (2010).
    The execution of model `fun` can be MPI parallelized.

    Parameters
    ----------
    fun : callable
        Function which computes the vector of residuals, with the signature
        ``fun(x, *args, **kwargs)``. The argument ``x`` passed to this
        function is an ndarray of shape (n,) (never a scalar, even for n=1).
        It must allocate and return a 1-D array_like of shape (n_obs,) or a scalar.
    x : array_like with shape (n,) or float
        Minimum point of independent variables obtained from POUNDerS or other routines. 
        If float, it will be treated as a 1-D array with one element.
    res_min: None or array_like with shape (n_obs,) or float (for n_obs==1), optional
        Vector of residuals at the minimum point. 
        If None (default), `fun` will be called to obtain this vector and `n_obs` will 
        be inferred from its length.
    bounds : 2-tuple of array_like, optional
        Lower and upper bounds on independent variables. Defaults to no bounds.
        Each array must match the size of `x` or be a scalar, in the latter
        case a bound will be the same for all variables. Use ``np.inf`` with
        an appropriate sign to disable bounds on all or some variables.
        For variables reaching the bounds, one-side finite difference method, instead
        of the two-side one, will be used for the derivative; they will also 
        be removed from the covariance / sensitivity matrix. 
    diff : array_like with shape (n,) or float
        Finite difference to be used to compute derivatives w.r.t. independent variables. 
        If float, all the variables will use the same difference. 
        For variables not reaching their bounds, derivatives w.r.t. them will be computed as
        (f(x+diff)-f(x-diff)) / (2*diff);
        otherwise, (f(x+diff)-f(x)) / diff or (f(x)-f(x-diff)) / diff to avoid exceeding 
        the bounds.  
    output_file: string or None, optional
        File path to write results. If None, no file output will be provided.  
    args, kwargs : tuple and dict, optional
        Additional arguments passed to `fun`. Both empty by default.
        The calling signature is ``fun(x, *args, **kwargs)``.
    comm: MPI communicator or None, optional
        If None (default), MPI_COMM_WORLD will be used. 
    
    Returns
    -------
    J: ndarray with shape (n, n_obs)
        Matrix of derivatives. The element of row i, column j is the derivative of the j-th residual
        w.r.t. the i-th variable. 
    cov: ndarray with shape (n, n)
        Covariance matrix. Variables reaching their bounds are removed. 
        None if all the variables reach their bounds. 
    eigvals: ndarray with shape (n,)
        Eigenvalues of the covariance matrix. 
        None if all the variables reach their bounds.
    eigvecs: ndarray with shape (n, n)
        Eigenvectors of the covariance matrix. 
        None if all the variables reach their bounds. 
    corr : ndarray with shape (n, n)
        Correlation matrix. Variables reaching their bounds are removed. 
        None if all the variables reach their bounds. 
    sen: ndarray with shape (n, n_obs)
        Sensitivity matrix. Variables reaching their bounds are removed. 
        None if all the variables reach their bounds.  
    
    References
    ----------
    .. M. Kortelainen, T. Lesinski, J. Moré, W. Nazarewicz, J. Sarich, N. Schunck, M. V. Stoitsov, and S. Wild,
       "Nuclear energy density optimization", Phys. Rev. C, 82, 024313 (2010).
    .. J. R. Donaldson and R. B. Schnabel, "Computational Experience With Confidence Regions and Confidence 
       Intervals for Nonlinear Least Squares", Technometrics 29, 67 (1987).
    """  
    
    # MPI Setup; different MPI processes have diffent values of "rank"
    if comm is None:
        if import_MPI:
            rank = MPI.COMM_WORLD.Get_rank()
        else: # sequential run
            rank = 0
    else:
        rank = comm.Get_rank()

    # prepare x
    assert(not np.iscomplexobj(x))
    x = np.atleast_1d(x).astype(float)
    assert(x.ndim == 1)
    n_param = len(x)

    # prepare diff
    diff = np.abs(np.asarray(diff, dtype=float))
    if diff.ndim == 0: diff = np.resize(diff, n_param)
    assert((diff.ndim == 1) and (len(diff) == n_param))

    # prepare bounds
    lb, ub = _unpack_tuple(bounds, n_param)
    assert(np.all(lb <= ub))
    reach_lb = (x<=lb)
    reach_ub = (x>=ub)
    reach_bound = reach_lb | reach_ub

    # prepare output file
    if rank == 0:
        try:
            outfile = open(output_file, 'w')
        except:
            printf('Warning from stat_analysis: no file output.')
            outfile = None
    else:
        outfile = None
    
    # prepare residual at minimum chi^2 
    if res_min is None:
        res_min = fun(x, *args, **kwargs)
    res_min = np.atleast_1d(res_min)
    n_obs = len(res_min)
    min_chi_sqr = np.sum(res_min**2)
    if outfile is not None:
        printf('# Parameter vector at min chi^2: \n', x.tolist(), file=outfile)
        printf('# Residual vector at min chi^2: \n', res_min.tolist(), file=outfile)
        printf('# Min chi^2: \n', min_chi_sqr, file=outfile)

    # finite difference method: go through all the varaibles
    derivatives = np.zeros((n_param, n_obs))
    for ind in reversed(range(n_param)): # reverse order, as special params are usually at the end
        flag = 0
        x_compute = deepcopy(x)
        if reach_lb[ind]:
            flag = 1
        if reach_ub[ind]:
            flag = -1
        if flag == 0:
            x_compute[ind] = x[ind] + diff[ind]
            res1 = fun(x_compute, *args, **kwargs)
            x_compute[ind] = x[ind] - diff[ind]
            res2 = fun(x_compute, *args, **kwargs)
            derivatives[ind, :] = (res1 - res2) / (2*diff[ind])
        else:
            x_compute[ind] = x[ind] + flag*diff[ind]
            res1 = fun(x_compute, *args, **kwargs)
            derivatives[ind, :] = (res1 - res_min) / (flag*diff[ind])
    if outfile is not None:
        printf('# Derivatives: \n', derivatives.tolist(), file=outfile)

    # statistical analysis
    if not np.all(reach_bound):
        J = derivatives[~reach_bound, :] # remove variables reaching bounds
        cov, eigvals, eigvecs, corr, sensitivity = stat_from_derivatives(min_chi_sqr, J, outfile)
    else:
        cov = eigvals = eigvecs = corr = sensitivity = None
        if outfile is not None:
            printf('# Warning: All the variables reach their bounds. No statistical analysis performed.', file=outfile)
    
    # finish
    if outfile is not None:
        outfile.close()
    return derivatives, cov, eigvals, eigvecs, corr, sensitivity  
    