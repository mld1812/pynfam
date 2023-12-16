# Overview

This repository contains the python package *pynfam* which includes the functions,
classes, and methods for running large-scale beta decay caluclations with the fortran
codes hfbtho and pnfam. It contains an executable script *run_pynfam.py* through which
the user can supply all necessary inputs for one or more beta decay calculations,
which are parallelized with mpi4py.

Currently, runs of fortran codes are handled by the *fortProcess* parent class in the
*fortran* sub-package. A run of an executable (i.e. the *runExe* method) consists of
creating a uniquely named directory, generating all necessary input files, then executing
the fortran program from within this directory using python's subprocess.Popen() method.
The fortran executables themselves are the only extra files needed to use this package.
(Future versions may replace this functionality with python bindings.)
A compatible copy of the hfbtho and pnfam source codes are included in this repository
for compiling their executables.

### Setup Quick Guide

1. Install the python package. If the setup.py file is located in the directory
   "/path/to/package", execute ``` pip install /path/to/package ```.
2. Compile the fortran code hfbtho using gnu make into executable "hfbtho_main".
3. Compile the fortran code pnfam using gnu make into executable "pnfam_main".
4. Move the executables into a common directory.
5. Use a copy of the script *run_pynfam.py* to execute serial or parallel calculations.
   Available inputs are listed in the pynfam/config.py file, and described in the
   pynfam_inputs and override_settings sections below.

### Execution Quick Guide

Each beta-decay calculation consists of:

1. HFBTHO runs for even-even, zero temperature HFB solutions with several deformations
2. * Feature: HFBTHO re-runs for non-converged even-even, zero temperature HFB solutions
3. HFBTHO runs for all blocking candidates of all even solutions, or for finite temperature
   HFB solutions.
4. Lowest energy HFB solution is chosen as the ground state and used as the basis for
   the pnFAM calculation.
5. pnFAM runs for each EQRPA point per operator.
6. Strength function is constructed from pnFAM results, phase space is calculated,
   and phase-space-weighted shape-factor is integrated for rates.

---

# Detailed Installation Guide

### Python package

Installation of the pynfam package can be performed with pip by calling install on the
directory which constains the setup.py file. E.g., if the path is "/path/to/pynfam/setup.py":

``` pip install /path/to/pynfam/ ```

If the user is not working in a virtual environemnt (see below) and does not have root
access to the python intsallation, doing ```pip install --user``` instead will install
the package locally, typically under ```~/.local/lib/pythonX.X/site-packages```.

### Dependencies

The *pynfam* package uses the following modules outside the standard library (all
can be found in the PYPI repository):

* numpy, pandas, scipy, mpi4py, f90nml, future

The package also depends on the HBFTHO and pnFAM executables. Compiling HFBTHO and
pnFAM requires the following libraries to be installed in the users environment.
One must change the appropriate Makefile variable to be the path to these libraries
when compiling the executables (see source for more details).

* LAPACK, BLAS

The openmp versions of the fortran programs should be compatible with the *pynfam* package
(though this has not been tested yet), but the versions compiled with MPI active should
NOT be used. All mpi parallelization is handled by python via mpi4py.

### mpi4py

A few errors encountered with mpi4py are worth noting here.

1. mpi4py must be installed and linked to the same mpi implementation currently active. This can
   cause issue e.g. on clusters where mpi4py was installed with one implementation but the user
   has a different one loaded in their modules. Sometimes the user must check their environment
   paths to ensure the correct MPI implementation will be found, and also mpi4py can be found by
   python (e.g. in the PYTHONPATH variable).
2. This package calls fork() in an mpi environment via python's subprocess package, which can be
   dangerous. Supposedly the versions of subprocess imported in this package are fork safe, and
   the fortran programs being executed do not themselves contain calls to MPI, so all should work
   fine. However, we have enountered and never were able to resolve some errors related to this
   using OpenMPI 3.0.0. Our recommendation is to use the latest version of mpi4py with mvapich2.

### Submodules

For pynfam users who have access to the primary HFBTHO repository, a git submodule to this
repository is included in the pynfam repository. A hard copy of the latest pynfam-compatible
HFBTHO repository version is included for users without this access.

After cloning the pynfam repository, the submodule directory will still be empty. To populate
the directory, users must run:
```git submodule init; git submodule update```
At this stage the the submodule will not be on any working branch (detached HEAD state), so
no changes will be tracked. To work on the project from the pynfam repository, be sure to
checkout a branch first.

For a short tutorial on using git submodules, see:
https://git-scm.com/book/en/v2/Git-Tools-Submodules

---

# Workflow

### Parallelization

Pynfam parallelizes within a single beta-decay calculations and among multiple beta-decay calculations.
The workflow for a single beta-decay calculation involves at most four serial stages that each
generate a list of "fortran objects" as tasks to be executed in parallel.

1. HFBTHO runs for even-even, zero temperature HFB solutions with several deformations
2. * Feature: HFBTHO re-runs for non-converged even-even, zero temperature HFB solutions
3. HFBTHO runs for all blocking candidates of all even solutions, or for finite temperature
   HFB solutions.
4. pnFAM runs for each EQRPA point per operator.

The MPI workflow is designed for load balancing across multiple beta-decay calculations, and offers
the best performance if the number of MPI processes is much larger than the number of beta-decay
calculations. MPI processes are split into two communicators, one with a small number of "master"
processes (<= number of beta-decay calculations), and another with a large number of "worker" processes.
A master process will execute pynfam code for a single beta-decay calculation, and multiple master
processes allows parallelization across multiple beta-decay calculations. As each master process
generates the task lists above, the tasks are distributed to the pool of worker processes, which is
shared by all master processes. In this way the work load is evenly distributed among the workers, even
if the number of tasks for different beta-decay calculations differs significantly.

This workflow requries one worker to be reserved as a "leader" to handle details of the
communications, therefore the resources are organized as:

* comm_size = nr_masters + nr_workers + 1_lead_worker

The code can also be run serially with a single process, but the MPI execution therefore requires
at least 3 processes. Some flexibility is provided for the user to determine the size of the
worker pool and the number of master processes via the input parameter nr_parallel_calcs
(see pynfam_inputs below). If the number of beta-decay calculations is greater than the number of
master processes, they are executed by the master processes in round-robin fashion.

### Runtime Behavior

There are several non-fatal error checks in pynfam which allow a large scale calculation
to continue while skipping an individual calculation which encountered the error. Examples
include a non-converged HFB solution, a segfault or other error encountered by the fortran
executable, etc.

---

# Inputs

Inputs for a beta-decay calculation are specified in the executable script *run_pynfam.py*.
A copy of the script is stored in the pynfam-generated output directory for future reference.
The contents of run_pynfam.py are split up into two main input sections, specified as the
dictionaries pynfam_inputs and override_settings which are discussed below.

### pynfam_inputs

The parameters in this section are common among all beta-decay calculations executed using
the same script.

* **directories (Str)**

    * 'outputs': Name of the directory for pynfam outputs
    * 'exes'   : Path to a single directory containing the HFBTHO and pnFAM executables
    * 'scratch': Path to location for pynfam outputs directory

* **nr_parallel_calcs (Int or None)**

     In the case of multiple beta-decay calculations (see override_settings for details on
     how to specify multiple calculations), each one may be computed in parallel. This
     parameter determines the number of "master" mpi processes that will be dedicated to running
     the set of beta-decay calculations. If less than the total number of calculations, the
     calculations will be executed by these processes in a round-robin fashion. The number
     of mpi worker processes dedicated to running tasks is comm_size - nr_parallel_calcs - 1.
     If nr_parallel_calcs is 0 or None, it defaults to the number of beta-decay calculations.

* **rerun_mode (Int)**

    This parameter only takes effect if the output directory already exists.

    * 0 - Full Restart: Uses the input file to determine the state of pre-existing data
          and pickup where the calculation left off if the data is incomplete. Non-converged
          HFB solutions count as incomplete, and will be re-run (unconstrained) to allow
          more iterations. Sub-directories must have original naming conventions.
    * 1 - FAM Restart: Determines the number of calculations based on number of existing directories.
          All existing directories must contain files for a completed hfb ground state, and the
          calculation begins in the FAM stage, ignoring any HFB inputs. Restart behavior for
          pre-existing FAM data is identical to Full Restart. Sub-directory names do not matter.
          Tuple override settings must have length equal to number of existing directories.
    * 2 - 2BC Restart: For calculating 2BC external field and storing the result. All HFB
          ground states must exist, this mode cannot be run from scratch. The parallelization
          scheme does not use any workers, only masters. This mode is meant to be run with
          openMP-enabled pnFAM, so each master can run a 2BC calculation on about 1 node of threads.
          In this way, operators are run serially, but nuclei are parallelized. See the README_tbc.md
          for more details.

* **gs_def_scan (Tuple)**

    This setting activates our built in method for finding axially deformed ground states.
    The input has the form (kickoff, (b2_1,b2_2,...)). For each value of deformation beta2
    (b2), hfbtho is run with the following inputs:

     * basis_deformation = b2
     * beta2_deformation = b2
     * lambda_active     = [0, kickoff, 0, 0, 0, 0, 0, 0]
     * lambda_active     = [0, Q2(A,b2), 0, 0, 0, 0, 0, 0]
     * Q2(A,b2)          = sqrt(5/pi)/100 * b2 * A^(5/3)

    The kickoff parameter maintains its original behavior as in hfbtho:

    * kickoff =  0 (or b2 iterable is empty): gs_def_scan feature is not used
    * kickoff = -1: Applies constraint for first few iterations only
    * kickoff = +1: Applies constraint for entire calculation

    Additionally it has options for blocking calculations of odd nuclei:

    * kickoff = +/-1: even calculations are constrained, odd are unconstrained
    * kickoff = -2: kickoff=-1 behavior for even and odd calculations
    * kickoff = +2: kickoff=+1 behavior for even and odd calculations

    Note: if active, these values will supercede inputs from 'override_settings'.

* **dripline_mode (Int)**

    Serially calculate out to the neutron dripline from the initial nucleus.

    * 0    - off
    * +/-1 - increase N by 1 until 1Sn = E(N+1)-E(N) < 0.
    * +/-2 - increase N by 2 until 2Sn = E(N+2)-E(N) < 0.

    The sign indicates direction to block. Blocking in the initial nucleus is
    independent of this sign.

* **ignore_nonconv (Int)**

    Define how to handle non-converged (nc) hfb solutions.

    * 0 - stop for any nc soln
    * 1 - discard odd nc solns (stops if ANY even or ALL odd are nc)
    * 2 - discard all nc solns (stops if ALL even or ALL odd are nc)
    * 3 - use nc solns

* **fam_contour (Str)**

    Specify the the type of energy contour to compute strength along. These come in 2 flavors:

    * Open contours: Used for computing strength functions. The imaginary part of the energy determines
      the half-width of the Lorentzian smeared strength peaks. Different open contours allow for a
      non-constant smearing, which gives greater flexibility in specifying the density of the energy
      grid and can save computational resources. Note that if the imaginary part is large,
      1/(x+iy) = P(1/x) + i /pi \delta(x) is no longer valid and the QRPA sum rule breaks.

    * Closed contours: Used for computing rates, but do not provide details of the strength function.
      This is done via a complex contour integration of the phase space weighted strength (or shape
      factor, which combines phase space weighted strength from allowed + first forbidden operators).

    All contours are defined on an energy interval along the real axis, while other settings specific
    to the type of contour define the behavior in the complex plane, spacing, number of points, etc.
    A contour's defining settings are listed in the config.py file, for details see the
    override_settings section below.

    A list of the available contours is given below. Contours with * adjust number of points until
    the desired profile fills the energy interval. Those without have a specified number of points.
    Where applicable, the spacing is defined in terms of the half-width via the de_hw_ratio setting.
    If de_hw_ratio > ~1, the contour grid is too coarse, and may miss peaks in the strength function.

    * None      : Skip the fam calculation, only compute the hfb ground state.
    * 'CIRCLE'  : (Closed) Circular contour
    * 'CONSTl'  : (Open)   Constant defined by nr_points (linspace)
    * 'CONSTr'  : (Open*)  Constant defined by spacing  (range)
    * 'EXP'     : (Open*)  Exponential profile
    * 'MONOMIAL': (Open*)  Monomial profile
    * 'FERMIs'  : (Open*)  Static Fermi function profile
    * 'FERMIa'  : (Open)   Adaptive Fermi function profile. Changes the contour itself to fill an
       energy interval with a specified number of points, prioritizing small width at small energy.
       This occurs in stages, between two limiting cases.

        1. Max points limiting case: CONSTl at a min half-width.
        2. Raise right side of contour to form a fermi function type profile.
        3. If right side reaches max half-width, raise left side of contour until max half-width.
        4. Min points limiting case: CONSTl at max half-width
        5. If the interval is still too large after stage 4, the number of points is increased until
           CONSTl at max half-width fills the interval.

* **beta_type (Str)**

    A string '-', '+', or 'c' for beta minus, beta plus, or electron capture. Note that the FAM strength
    for '+' and 'c' are identical, but the beta solution will differ.

* **fam_ops (Str or Tuple)**

    pnfam calculates the hfb linear response to allowed and first-forbidden beta decay operators
    (external fields). The operators are separated by how they change K=(-J,...,J), thus each
    (Operator, K) combination will have it's own strength function. For even-even parents, and odd
    parents computed in the equal filling approximation, only K>=0 contributions need to be computed.

    The following operators are implemented in pnfam:

    * 'F'   : (J=0) Fermi operator (identity)
    * 'GT'  : (J=1) Gamow-Teller operator (pauli matrix sigma)
    * 'R'   : (J=1) position vector (r)
    * 'P'   : (J=1) momentum over i (p/i)
    * 'RS0' : (J=0) r and sigma coupled to J=0
    * 'RS1' : (J=1) r and sigma coupled to J=1
    * 'RS2' : (J=2) r and sigma coupled to J=2
    * 'PS0' : (J=0) p/i dotted with sigma (note: this is not equal to p/i and sigma coupled to J=0)

    Each operator is accompanied by iso-spin raising or lowering operator, dicated by the beta_type input.
    This package wraps sets of operators for convenient beta decay calculations. Available calculations
    include total, total allowed, total first-forbidden, or transition by J-pi. Calculation names,
    number of strength functions computed (n), and operators considered, are listed below.
    (note only K>=0 operators are computed):

    * 'All'       : n=14, all operators
    * 'Allowed'   : n= 3, [F,GT]
    * 'Forbidden' : n=11, [R,P,RS0,RS1,RS2,PS0]
    * '0+'        : n= 1, [F]
    * '1+'        : n= 2, [GT]
    * '0-'        : n= 2, [RS0,PS0]
    * '1-'        : n= 6, [R,P,RS1]
    * '2-'        : n= 3, [RS2]
    * (Operator,K): n= 1, manually specify 1 strength function with a tuple of form (str, int)
      where str is one of the 8 operators listed above and int is the K value.

### override_settings

Default calculation settings are defined in the config.py file in the pynfam package, and
can be overridden here. To run multiple calculations in parallel, a tuple of values can be
specified for any parameter. This will run len(tuple) beta-decay calculations, with the i'th
calculation using the i'th tuple element as the setting. Standard (Non-tuple) parameters are
applied to every calculation. All tuple inputs must have the same length.

* **hfb (Dict)**

    HFBTHO namelist parameters to override.

* **ctr (Dict)**

    Parameters defining the fam_contour. All contours have common settings `energy_min`,
    `energy_max`, and `hfb_emin_buff`. The latter setting defines a buffer to increase the default
    interval computed from properties of the HFB ground state by reducing the default value of
    energy_min (see below). If energy_min is specified, this value will supercede the above behavior.

    The default energy interval computed from the HFB ground state is:

    * energy_min = min(0, E_gs - hfb_emin_buff)
    * energy_max = (E_gs + Q)

        - E_gs_EE = min(E_p+E_n)
        - E_gs_EO = min(E_p - E_n_blocked), min(E_n-E_p_blocked)
        - E_gs_OO = - En_blocked - Ep_blocked
        - Q_bm ~ MnH - lambda_p + lambda_n - E_gs
        - MnH = M_neutron - M_Hydrogen
        - lambda = HFB approximation to the Fermi energy

    If an imaginary value is supplied as the `half_width` setting for applicable contours, strength
    is on  [energy_min*1j, energy_max*1j] a distance Imag(half_width) from the imaginary axis.

* **fam (Dict)**

    pnFAM namelist parameters to override. The default values include EDF coupling overrides
    for values fit by Mustonen 2016 (Table 2 set 1A from PRC 93 014304). To use the pnfam code's
    default values (computed from the interaction) the user can override these values by
    supplying the empty string ''. This leaves the namelist entry for the parameter blank.

* **psi (Dict)**

    Parameters related to the phase space and shape factor. For every QRPA energy, there is a
    phase space integral that is calculated with gauss-legendre quadrature. The integrals as a
    function of QRPA energy must be approximated by an analytic function for closed contours.
    There are two main methods for approximating the phase space integrals in the complex plain:
    rational function interpolation (RATINT) and polynomial least squares fit (POLYFIT).

    Since FAM and HFB do not depend on phase-space related parameters, and psi-related calculations are not
    computationally expensive or parallelized, we allow passing a list for any psi setting. Lists act similarly
    to tuples, but are computed WITHIN the beta_soln directory for a single beta-decay calculation. The FAM
    contour can sometimes depend on the settings Q_eff and log(pYe). In these special cases, the different
    contours are also computed WITHIN the fam_soln directory.

    * psi_approx (Str): Either RATINT or POLYFIT to determine which set of settings is used to approximate
      an analytic function of phase space integrals. Ignored for open contours.
    * psi_glpts (Int): Number of gauss-legendre nodes used to compute the phase space integrals.
    * screening (Bool): Apply Rose screening to account for screening from atomic electrons. Only applied
      for beta-plus decay.
    * Q_eff (Float): Effective Q-value to use for phase-space. Changes the QRPA energy where phase space is
      zero from the default EQRPA_max=Q_hfb+E_gs.
    * Q_eff_mode (Int): Takes values +/-[0,1,2]. 0 means Q_eff is not used. 1 uses Q=Q_eff.
      2 uses Q=Q_hfb+Q_eff. Minus sign means phase space is shifted, but every Q_eff uses
      same contour with default EQRPA_max. Positive sign means both phase space and EQRPA_max are shifted.
      This means for closed contours 1 contour per Q_eff is calculated, while for open contours 1 contour
      is calculated out to the largest EQRPA_max determined from any Q_eff.
    * GA (Float): Value of the axial vector coupling used in the shape factor.
    * GV (Float): Value of the vector coupling used in the shape factor.
    * log(pYe) (Float): Log of electron fraction times stellar density. Only applies to
      finite temperature electron capture calculations. List computes rates for multiple densities
      using the same FT-QRPA strength function. For closed contours, supplying a list might use a
      different contour, or two contours, for a given density. This is for numerical stability of the
      analytic phase space approximation.
    * ratint_pts (Int): Number of points (on a Chebychev grid) used for the rational function interpolation
      of phase space integrals. Only used for closed contours with psi_approx = RATINT.
    * polynomial_order (Int): Order of the fitting polynomial for phase space integrals. Only used for
      closed contours with psi_approx = POLYFIT.
    * polynomial_fit_pts (Int): Number of points (on a Chebychev grid) used for the polynomial fit for phase
      space integrals. Only used for closed contours with psi_approx = POLYFIT.

---

# Finite Temperature

### The prefactor

The pnFAM response function at finite temperature contains ensemble averaged strength of the forward
and reverse processes (the latter comes with negative sign). Moreover, it includes transitions from
all excited parents to all excited daughters, which includes negative energy transitions (de-excitations).
Because of this, all forward and reverse strength overlaps. We must used detailed balance to remove
the reverse process strength, which amounts to introducing a prefactor S_FT(w)=(1/(1-exp(-w/KT))S_FAM(w).
This prefactor increases strength at w>0, account for the overlap of negative reverse process strength,
and also makes strength at w<0 positive by the correct amount. However, it contains a pole exactly at w=0.

Python handles the prefactor, not pnfam. This is because we pretend that strength on open contours is
on the real axis, and thus w in the prefactor can be treated as real. Closed contours need complex omega.
Since pnfam2 has no concept of the contour, this is handled when combining all strength points and writing
the fam_soln results. So, fam_soln results include the prefactor, but fam_meta results do not!

We must avoid the pole at zero from the FT prefactor for closed contours. Technically, there are many poles
along the w=0 axis, at 2pi*n*i. So it's best for contour calculations to be split up into one contour at w<0
and one at w>0, avoiding the w=0 axis altogether. Thus if the specified energy range includes zero, two contours
will be computed, in addition to a third contour around just the pole. It turns out the contours very close to
the pole do not give the correct result, unless the points are rotated in a very particular way that I don't yet
understand. However, if the pole is actually ON the contour, as long as no points on the contour are
exactly on top of the pole, we very reliably get an integral that contains 1/2 the residue of the pole
(one can prove this with calculus of paritial residues). So the approach is the compute three contours,
[-E, 0]+[0,+E]+[-eps,eps], where the first two actually touch zero and each pickup 1/2 the residue, and the third
picks up the full residue of just the pole, which can then be subtracted out.

IMPORTANTLY, results written to beta_soln already have 1/2 the pole residue subtracted out! This makes it
more clear the meaning of zeroed_beta.out, since the pole calculation is not expected to be positive.

### Energy Range

FT calculations need strength from [-infty,EQRPAmax] for spontaneous decay, and [-inft,+infty]
for FT electron capture, since energetic electron can still be captured by nuclei with negative Q values.
The lower bound can be determined by the FT prefactor, which rapidly dies to zero at very negative energies.
The upper bound for electron capture can be determined by the phase space, which dies off rapidly at large energies. 
This is how default bounds are determined, up to an arbitrary maximum of [-30,+30]MeV.

Note that circles extending far into the complex plane behavior poorly, particularly due to the oscillations in the
phase space approximations, so crank up the interval with care. Also, convergence of pnfam at high energies becomes poor.
An ellipse might allow computing to higher energies without going far into the compelx plane, but be careful because
the pole residue is only accurately computed using a small circle.

Finally, EQRPAmax should be determined using the zero-temperature fermi energies. Thus, for each FT calcualtion
we must also compute the 0T hfb solution, which is included as an additional subdirectory.

### FT Electron Capture Phase Space

The FT electron capture phase space requires two numerical integrations, a Gauss-Legendre to capture an initial
hump near the electron fermi energy, and a Gauss-Laguerre to integrate an exponentially decaying tail to infinity.
It typically requires larger number of psi_glpts than other modes, around 60-100, for best accuracy.

Also, the FT electron capture phase space changes from mostly exponential to mostly polynomial at the
electron fermi energy (EQRPA(mu) = mu + EQRPAmax). For the analytic approximation to do well, we
must split the contour into two if the fermi energy falls inside the contour, using a either an exponential
or a polynomial fit. RATINT will not fit the function well and WILL introduce spurious poles! Note that the
fermi energy is a function of temperature and log(pYe), so the split may be different depending on densities.
This means, at most, for 1 density we may need 3 contours, pluts the pole contour, for 4 total. The code
does it's best to identify non-unique contours and avoid computing them, which is handy when requesting
multiple densities at once.

The phase space dies off rapidly to zero, and the contour integral and/or phase space approximation will become
inaccurate if a large region of zero is included. Thus the contour right bound is adjusted to something less than
or equal to the request value to avoid this.

---

# Outputs

The output directory tree is structured as detailed below, where parenthesis indicate files.
pynfam generates logfiles along the way with key outputs. Key hfb+fam+beta outputs for every
beta-decay calculation are compiled into one master logfile in the meta data directory.

* For calculations with lists of psi inputs the contents of beta_soln and beta_meta directories
  contain one sub-directory per element in the lists, labelled psi_xxx, where xxx is the element's
  index in the list. Note that the subdirectory may be empty if the phase space is exactly zero
  at all QRPA energies.

* For FT calculations where a closed contour is split, additional subdirectories will be created.

* For lists of log(pYe) or Q_eff (with Q_mode > 0), the fam_soln and fam_meta directories will contain
  one sub-directory per unique contour calculated.

        outputs
          |- meta
          |    |- (copied_inputs_and_slurm_files+master_logfile)
          |
          |- calc_label
               |- hfb_soln
               |    |- (copied_ground_state_output_files+mini_logfile)
               |    |- hfb_meta
               |         |- (logfile_for_all_hfb_runs)
               |         |- (tarfile_of_all_hfb_run_dirs)
               |              |- hfb_label_even
               |              |    |- (raw_output_files)
               |              |- hfb_label_odd
               |                   |- core
               |                       |- (raw_output_files)
               |                   |- candidate_label
               |                       |- (raw_output_files)
               |
               |- fam_soln
               |    |- (strength_output_text_and_binary_files_per_operator)
               |    |- fam_meta
               |         |- (representative_logfile_for_pnfam_runs)
               |         |- (tarfiles_of_all_pnfam_runs_per_operator)
               |              |- fam_point_label
               |                   |- (raw_output_files)
               |
               |- beta_soln
                    |- (beta_decay_output_files+mini_logfile)
                    |- beta_meta
                         |- (non_essential_output_files)
