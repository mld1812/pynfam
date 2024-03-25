r"""
    Author: Tong Li, Mich State Univ, 2021-
"""
r"""
    Notes on the installation of petsc and petsc4py:
    Run command "python -m pip install petsc petsc4py" to install.
    Add option "--user" if no administration privilege is granted.
    BLAS and LAPACK should be present (loaded) for the installation.
    The module OpenBLAS / Intel MKL contains both BLAS and LAPACK.
    The default modules on MSU iCER should work well. 
"""

import numpy as np
import sys
from copy import deepcopy
import functools
printf = functools.partial(print, flush=True) # force output flush
try:
    # expose stat_analysis to scripts that import pounders
    from stat_chi2opt import stat_analysis, stat_from_derivatives
except ImportError:
    printf('Warning: Cannot import stat_analysis or stat_from_derivatives! Fall back to empty functions.')
    stat_analysis = lambda *args: None
    stat_from_derivatives = lambda *args: None
try:
    from mpi4py import MPI
    import_MPI = True
except ImportError:
    import_MPI = False

class _comm_seq(object):
    # MPI communicator for the sequential case
    def Get_rank(self): 
        return 0
    def bcast(self, data, root=0): 
        return data

class _termination_flag(object):
    # termination flag, for MPI model execution
    pass

_reason_map = { 3: "||G(X)|| < gatol", \
                4: "||G(X)|| / |F(X)| < grtol", \
                5: "||G(X)|| / ||G(X0)|| < gttol", \
                6: "Step size small", \
                7: "F < F_min", \
                8: "User defined convergence", \
                -2: "Iteration number limit exceeded", \
                -4: "Numerical problems", \
                -5: "Function evaluation number limit exceeded", \
                -6: "Line search failure", \
                -7: "Trust region failure", \
                -8: "User defined divergence"}

def _unpack_tuple(tin, n_param):
    assert(len(tin)==2)
    lb, ub = [np.asarray(b, dtype=float) for b in tin]
    if lb.ndim == 0: lb = np.resize(lb, n_param)
    if ub.ndim == 0: ub = np.resize(ub, n_param)
    assert((lb.ndim == 1) and (ub.ndim == 1) and (len(lb) == len(ub) == n_param))
    return lb, ub


def pounders_seq_wrapper(fun, x0, n_obs=None, bounds=(-np.inf,np.inf), scale_intervals=(None,None),
                         gatol=None, grtol=None, gttol=None, max_iter=None, max_nfev=None,
                         pounders_args=None, history_file='pounders.txt', args=(), kwargs={}):
    r"""Solve a nonlinear least-squares problem with bounds on the variables 
    using POUNDerS (serial, without MPI parallelization) provided by PETSc/TAO,
    with an interface similar to scipy.optimize.least_squares.

    Given the residuals f(x) (an m-D real function of n real
    variables), `pounders_seq_wrapper`
    finds a local minimum of the cost function F(x)::
        minimize F(x) = sum(f_i(x)**2), i = 0, ..., n_obs - 1),
        subject to lb <= x <= ub.

    A similar wrapper function is provided in Python package "Estimagic",
    a package for economic models.

    Parameters
    ----------
    fun : callable
        Function which computes the vector of residuals, with the signature
        ``fun(x, *args, **kwargs)``, i.e., the minimization proceeds with
        respect to its first argument. The argument ``x`` passed to this
        function is an ndarray of shape (n,) (never a scalar, even for n=1).
        It must allocate and return a 1-D array_like of shape (n_obs,) or a scalar.
    x0 : array_like with shape (n,) or float
        Initial guess on independent variables. If float, it will be treated
        as a 1-D array with one element.
    n_obs: int or None, optional
        Number of observations, or the dimension of the output of `fun`. 
        If None (default), `fun` will be evaluated once to obtain the dimension,
        which is not recommended if the evaluation is time-consuming.
    bounds : 2-tuple of array_like, optional
        Lower and upper bounds on independent variables. Defaults to no bounds.
        Each array must match the size of `x0` or be a scalar, in the latter
        case a bound will be the same for all variables. Use ``np.inf`` with
        an appropriate sign to disable bounds on all or some variables.
    scale_intervals : 2-tuple of array_like, optional
        Lower and upper bounds of scaling intervals. Defaults to no scaling.
        Variables will be linearly mapped from their scaling intervals to a unit  
        hypercube of [0,1]^n. Each array must match the size of `x0` or be a
        scalar, in the latter case a bound will be the same for all variables. 
        Use ``np.nan`` or None to disable scaling intervals on all or some variables.
    gatol: float or None, optional
        One of convergence criteria for POUNDerS: ||G(X)|| < gatol, 
        where G is the gradient of F. 
        If None (default), the default value in PETSc/TAO (1e-8) will be used. 
        Set it as zero to disable this option. 
    grtol: float or None, optional
        One of convergence criteria for POUNDerS: ||G(X)|| / |F(X)| < grtol, 
        where G is the gradient of F. 
        If None (defulat), the default value in PETSc/TAO (1e-8) will be used. 
        Set it as zero to disable this option. 
    gttol: float or None, optional
        One of convergence criteria for POUNDerS: ||G(X)|| / ||G(X0)|| < gttol, 
        where G is the gradient of F and G(X0) is the initial gradient.
        If None (defulat), the default value in PETSc/TAO (0.0) will be used. 
        Set it as zero to disable this option. 
    max_iter : int or None, optional
        Maximum number of POUNDerS iterations before the termination.
        If None (default) no limit is set; 
        but POUNDerS has a large upper limit (2000) hard coded for the maximum 
        number of iterations, which can be overwritten by `pounders_args`.
    max_nfev: int or None, optional
        Maximum number of function evaluations before the termination.
        If None (default) no limit is set; 
        but POUNDerS has a large upper limit (4000) hard coded for the maximum 
        number of function evaluations, which can be overwritten by `pounders_args`.
    pounders_args: string, tuple or list of strings, or None, optional
        Advanced settings for POUNDerS given by command line arguments.
        If None (default) no advanced setting will be used. Read TAO Users Manual for details.
    history_file: string or None, optional
        File path to write the history of POUNDerS iterations. 
        If None, no output file will be produced. 
    args, kwargs : tuple and dict, optional
        Additional arguments passed to `fun`. Both empty by default.
        The calling signature is ``fun(x, *args, **kwargs)``.

    Returns
    -------
    x : ndarray, shape (n,)
        Solution x for the minimization of F(x).
    summary: dict
        A summary of the solution status with following fields.
        * iter: int.
            Iteration number.
        * sol_scaled: ndarray, shape (n,)
            Solution x scaled w.r.t. `scale_intervals`.
        * residual: ndarray, shape (n_obs,)
            Vector of residuals at point x.
        * Fval: float.
            Value of F(x). 
        * gnorm: float. 
            Norm of the gradient of F(x).
        * reason: int.
            Iteration stopping reason. See PETSC/TAO manual for details.
        * reason_message: string. 
            Iteration stopping reason given in a human-readable format. 
    
    Attention
    ---------
    Default values used in POUNDerS may be different in different versions of PETSc,
    and might be changed in the future. 
    The values given above are taken from the source code of PETSc 3.15.1. 

    References
    ----------
    .. A similar implementation in Python package "estimagic"
       https://estimagic.readthedocs.io/en/v0.0.29/_modules/estimagic/optimization/pounders.html
    .. PETSc/TAO Web page and Users Manual
       https://petsc.org
    .. M. Kortelainen, T. Lesinski, J. Moré, W. Nazarewicz, J. Sarich, N. Schunck, M. V. Stoitsov, and S. Wild,
       "Nuclear energy density optimization", Phys. Rev. C, 82, 024313 (2010).
    .. S. M. Wild, J. Sarich and N. Schunk, "Derivative-free optimization for parameter estimation 
       in computational nuclear physics", J. Phys. G: Nucl. Part. Phys., 42, 034031 (2015).
    .. S. M. Wild, "Solving derivative-free nonlinear least squares problems with POUNDERS",
       Advances and Trends in Optimization with Engineering Applications (SIAM, Philadelphia, 2017),
       pp. 529–540.
    """

    # initialize PETSc with settings given by pounders_arg
    assert((pounders_args is None) or (isinstance(pounders_args, (tuple,list,str))))
    import petsc4py
    if import_MPI:
        comm_self = MPI.COMM_SELF
    else:
        comm_self = None
    petsc4py.init(args=pounders_args, comm=comm_self)
    from petsc4py import PETSc

    # prepare x0
    assert(not np.iscomplexobj(x0))
    x0 = np.atleast_1d(x0).astype(float)
    assert(x0.ndim == 1)
    n_param = len(x0)

    # prepare bounds
    lb, ub = _unpack_tuple(bounds, n_param)
    assert(np.all((lb <= ub) & (x0 >= lb) & (x0 <= ub)))
    lb = PETSc.Vec().createWithArray(np.where(np.isneginf(lb), PETSc.NINFINITY, lb), comm=comm_self)
    ub = PETSc.Vec().createWithArray(np.where(np.isposinf(ub), PETSc.INFINITY, ub), comm=comm_self)

    # prepare scaling
    ls, us = _unpack_tuple(scale_intervals, n_param)
    ds = us - ls
    def scale_to_unit_cube(xin):
        return np.where(np.isnan(ls) | np.isnan(us), xin, (xin-ls)/ds)
    def scale_from_unit_cube(xin):
        return np.where(np.isnan(ls) | np.isnan(us), xin, ls+ds*xin)

    # check max_iter and max_nfev, and prepare n_obs
    assert(((max_iter is None) or (max_iter > 0)) and ((max_nfev is None) or (max_nfev > 0)))
    if n_obs is None:
        res = np.atleast_1d(fun(x0, *args, **kwargs))
        n_obs = len(res)

    # prepare PETSc vectors
    x = PETSc.Vec().createWithArray(scale_to_unit_cube(x0), comm=comm_self)
    f = PETSc.Vec().create(comm_self)
    f.setSizes(n_obs)
    x.setFromOptions()
    f.setFromOptions()

    # residual function wrapper
    res_min = [np.inf, None]
    nfev = 0
    def formResidual(tao, xin, fout):
        nonlocal nfev, res_min
        temp_array = np.atleast_1d(fun(scale_from_unit_cube(xin.array), *args, **kwargs))
        temp_array = np.where(np.isneginf(temp_array), PETSc.NINFINITY, temp_array) # deal with -inf
        fout.array = np.where(np.isposinf(temp_array), PETSc.INFINITY, temp_array) # deal with +inf
        nfev += 1
        chi_square = np.sum(fout.array**2)
        if chi_square < res_min[0]:
            res_min[0] = chi_square
            res_min[1] = deepcopy(fout.array)

    # convergence test for max_iter and max_nfev
    def test_max_iter_nfev(tao):
        if max_iter is not None: 
            if tao.getIterationNumber() >= max_iter:
                return -2
        if max_nfev is not None:
            if nfev >= max_nfev:
                return -5
        return tao.getConvergedReason()

    # initialize POUNDerS
    tao = PETSc.TAO().create(comm_self)
    tao.setType(PETSc.TAO.Type.POUNDERS)
    try:
        tao.setResidual(formResidual, f)
    except AttributeError:
        # for compatibility with old version
        tao.setSeparableObjective(formResidual, f)
    tao.setTolerances(gatol=gatol, grtol=grtol, gttol=gttol)
    tao.setVariableBounds((lb, ub))
    tao.setConvergenceTest(test_max_iter_nfev)
    tao.setFromOptions()
    
    # set a monitor
    def monitor_func(tao, history, prefix):
        x = tao.getSolution().array
        its, Fval, gnorm, _, _, _ = tao.getSolutionStatus()
        printf(prefix+str(its), x.tolist(), scale_from_unit_cube(x).tolist(), \
              Fval, gnorm, file=history, sep='; ')
    try:
        history = open(history_file, 'w')
        printf('iter    x_scaled    x    chi^2     ||g||', file=history)
        prefix=''
    except:
        printf('Warning from pounders_seq_wrapper: history file cannot be opened or written. Fall back to stdout.')
        history = sys.stdout
        prefix='POUNDerS '
    tao.setMonitor(monitor_func, (history, prefix))

    # solve
    tao.solve(x)

    # result
    sol_scaled = deepcopy(tao.getSolution().array)
    sol = scale_from_unit_cube(sol_scaled)
    its, Fval, gnorm, _, _, reason = tao.getSolutionStatus()
    if reason in _reason_map:
        reason_message = _reason_map[reason]
    else:
        reason_message = "Undefined reason"
    summary = {'iter':its, 'sol_scaled':sol_scaled, 'residual':res_min[1], 'Fval':Fval, 'gnorm':gnorm, \
               'reason':reason, 'reason_message': reason_message}

    # clean
    if history_file is not None: history.close()
    lb.destroy(); ub.destroy()
    x.destroy(); f.destroy()
    tao.destroy()

    return sol, summary


def pounders_mpi_wrapper(fun, x0, n_obs=None, bounds=(-np.inf,np.inf), scale_intervals=(None,None), 
                         gatol=None, grtol=None, gttol=None, max_iter=None, max_nfev=None, 
                         pounders_args=None, history_file='pounders.txt', args=(), kwargs={}, comm=None):
    r"""Solve a nonlinear least-squares problem with bounds on the variables 
    using POUNDerS provided by PETSc/TAO.
    The execution of model `fun` can be MPI parallelized.
    See the documentation of `pounders_seq_wrapper` for details of parameters and returns.
    The only new parameter is `comm`, the MPI communicator;
    if None (default), MPI_COMM_WORLD will be used. 
    """  

    # MPI Setup; different MPI processes have diffent values of "rank"
    if comm is None:
        if import_MPI:
            comm = MPI.COMM_WORLD
        else: # sequantial run 
            comm = _comm_seq()
    rank = comm.Get_rank()

    # wrapper for fun
    def fun_wrapper(xin, *args, **kwargs):
        x = comm.bcast(xin, root=0)
        if isinstance(x, _termination_flag): # receive termination flag
            return x
        return fun(x, *args, **kwargs)

    # main work
    if rank == 0:
        tflag = _termination_flag()
        r = pounders_seq_wrapper(fun_wrapper, x0, n_obs, bounds, scale_intervals, \
                                 gatol, grtol, gttol, max_iter, max_nfev, pounders_args, \
                                 history_file, args, kwargs)
        comm.bcast(tflag, root=0) # broadcast termination flag
    else:
        while True:
            r = fun_wrapper(None, *args, **kwargs)
            if isinstance(r, _termination_flag): # receive termination flag
                break
    r = comm.bcast(r, root=0) # broadcast final result
    return r
