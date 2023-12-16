The pnFAM Quick Reference
=========================

There are currently two pnFAM executables in the project:  The one named
*pnfam_nompi.x* can be compiled to run in serial or with OpenMP, and the one
named *pnfam_mpi.x* can be compiled as purely MPI or as a hybrid
MPI+OpenMP code.  The latter can use ADLB for load-balancing, and it can
process several input files in parallel.
 
 
Installation
------------

### Dependencies

For input, the pnFAM codes require a binary output file from a modified version of
HFBTHO v2.00d. The required modifications to the version published in the CPC archive
are included with the pnFAM as a patch file.
Both the MPI and non-MPI versions need to be linked to the Gnu Scientific Library (GSL).
The MPI version additionally requires the ADLB load-balancing library, available for free
from https://www.cs.mtsu.edu/~rbutler/adlb/ .
The programs have been tested to successfully compile with sufficiently new gfortran
and Intel Fortran compilers on various platforms ranging from an Apple MacBook Pro
running OS X Yosemite to the Stampede supercomputer in Texas Advanced Computing Center
(TACC) and Hopper in National Energy Research Scientific Computing Center (NERSC).

### Step-by-step installation instructions

1. Download HFBTHO version 2.00d, available in the Computer Physics Communications
   Program Library (http://cpc.cs.qub.ac.uk).
2. Copy the file hfbtho.patch to the HFBTHO directory and then apply our patch that
   extends HFBTHO output using the command `patch < hfbtho.patch`
   in the HFBTHO directory.
3. Modify the HFBTHO Makefile as needed for your platform and compile HFBTHO.
4. (Needed only for compiling the MPI version:) Download and compile ADLB
   (https://www.cs.mtsu.edu/~rbutler/adlb/) in a separate directory, following their
   instructions. Then either copy the files _adlbf.h_, _libadlb.a_, and _libfadlb.a_ over to
   the pnFAM directory, or modify the pnFAM makefile so that they are found.
5. Modify the pnFAM makefile for the compiler installed on your platform, if needed, and
   compile.
6. (Recommended:) Run the Bash script `test_nompi.sh` in the tests/ subdirectory to confirm
   the non-MPI codes were properly compiled. You will likely have to set the paths in
   the beginning of the script to reflect your file organization for the tests to pass.


Running in a single-file mode
-----------------------------

In a single-file mode, the parameters are read by default from file *default.in*.
The input parameters are interpreted as follows:

### Namelist *GENERAL*

- *fam_mode*: currently the options are "STR" for computing strength functions,
"CONTOUR" for computing closed contours for beta-decay rates, and "FINDMAX" for
locating a (local) maximum of the strength function (useful for finding
resonance peak energies).  No default, must be given explicitly.
- *hfb_input_filename*: the file containing the special binary output from the
HFBTHO designed for the pnFAM.
- *fam_output_filename*: the final output file for pnFAM.
- *log_filename*: a file where verbose information about the run should be
printed.  If empty, writes to standard output instead.

### Namelist *EXT_FIELD*

- *operator_name*: selected operator to compute.  Currently implemented options
  are:
    + 'F+', 'F-': Fermi operator for beta+ and beta-
    + 'GT+', 'GT-': Gamow-Teller operator for beta+ and beta-
    + 'RS0+', 'RS0-': r and sigma coupled to J=0
    + 'PS0+', 'PS0-': p/i-dot-sigma (note: this is not equal to p/i and sigma coupled to J=0)
    + 'RS1+', 'RS1-': r and sigma coupled to J=1
    + 'R+', 'R-': r for beta+ and beta-
    + 'P+', 'P-': p/i for beta+ and beta-
    + 'RS2+', 'RS2-': r and sigma coupled to J=2
    
  Please note that the computation of the decay rate currently only supports
  beta- operators, while the strength functions can be computed for both.

- *K*: the angular momentum projection to the intrinsic axis.  Since parity
  is a good quantum number in the basis, negative values yield the same result
  as positive ones.
- *compute_crossterms*: if .true., also computes the cross-terms of the selected
  operator needed to evaluate non-unique forbidden decays

### Namelist *INTERACTION*

- *interaction_name*: selection of a built-in interaction or an external
  coupling constant file.  An external file 'myinteraction.dat' would be selected
  by specifying 'FILE:myinteraction.dat' as *interaction_name*.  The residual
  interaction can be turned off completely by selecting 'NONE'.
- *force_j2_terms*: if .true., the J^2 terms are included in the computation,
  even if they shouldn't be ordinarily (e.g. forces such as SLy4 or SkM* where
  the J^2 terms are neglected in the HFB computation).  If left .false., the
  J^2 terms are chosen to match what was used in the HFB computation (e.g. for
  SLy4 there are no J^2 terms, but for SLy5 they are included).
- *require_gauge_invariance*: if .true., the program (1) adjusts the time-odd
  coupling constants CT and CF to obey the gauge relations linking them to the
  J^2 couplings and then (2) aborts if any coupling constants do not obey their
  individual gauge relations.  Used to guard against unintentional breaking of
  the gauge invariance.
- *require_self_consistency*: if .true., the program aborts if the coupling
  constants of the residual interaction are not consistent with those of the
  HFB solution read from the file.  Used to guard against unintentional
  breaking of self consistency of the interaction.
- *vpair_t0*: overrides the default isoscalar pairing strength
- *vpair_t1*: overrides the default isovector pairing strength
- *override_cs0*, *override_csr*: overrides the default density-dependent
  coupling of the spin density term of the EDF
- *override_cds*: same for the laplacian s term
- *override_ct*: same for the T term
- *override_cf*: same for the F term
- *override_cj*: same for the j term
- *override_cgs*: same for the div s term
- *override_csdj*: same for the curl j term

### Namelist *STR_PARAMETERS*

This namelist is only read if *fam_mode* is set to "STR".

- *energy_start*: the starting value for the energy
- *energy_step*: the spacing of energy values
- *nr_points*: the number of energy values to compute.
- *half_width*: the imaginary part added to the
  energies to introduce smearing of the strength function.  Note that smaller
  half-width requires more closely spaced points to be computed (i.e. smaller
  *energy_step*) and that in general, the solution of the non-linear pnFAM
  equations takes more iterations when closer to the real axis (where the poles
  of the strength function lie)

### Namelist *FINDMAX_PARAMETERS*

This namelist is only read if *fam_mode* is set to "FINDMAX".  In the FINDMAX
mode, the code first moves "uphill" from the initial guess until it finds an
interval where the sign of the derivative flips, and then bisects until the
interval where the sign flip occurs is as small as requested.

- *energy_start*: the initial guess for the local maximum.
- *search_step*: the initial step size for locating the interval where the maximum
  is (where the derivative flips sign).
- *max_points*: the maximum number of points to try before giving up.  Note
  that since two points are evaluated for approximating the derivative at each
  point, the actual maximum number of strength function evaluations is twice this.
- *half_width*: the half-width for the strength function.
- *delta*: the points computed for approximating the derivative at energy _x_ 
  are _x+delta_ and _x-delta_.
- *energy_tolerance*: the requested accuracy for the location of the maximum.

### Namelist *CONTOUR_PARAMETERS*

This namelist is only read if *fam_mode* is set to "CONTOUR".

- *q_value*: the Q value of the decay.  The
  contour in the decay rate code is a circle that crosses the real axis at
  the origin and *q_value* + ground state energy; the ground state energy is
  approximated with the sum of the lowest neutron and proton quasiparticle
  energies. In the absence of pairing correlations, the ground-state energy is
  the lowest neutron hole-state energy plus the lowest proton particle-state
  energy.
- *nr_points*: the number of energy values to compute.  The
  actual number of computed points is about half of this because symmetry is
  used to reduce the number of needed pnFAM computations.
- *rot_correction*: rotational energy correction to
  approximate the mismatch between the intrinsic and lab systems, in MeV.
  If negative, computed from the HFB solution using the cranking formula adapted
  for the HFB.  If zero, neglected.
- *include_endpoint*: if *.true.*, explicitly computes the strength function
  in the endpoint where it usually does not contribute due to the vanishing
  phase-space.
- *rotate_contour*: the continuous circular contour we use is discretized into
  *N* points and *N-1* intervals of angle *t_N*.  When *.true.*, this option
  avoids driving through a pole on the *Re(omega)* axis by rotating each of
  the *N* contour points by an angle *t_N / 2*.  This is helpful primarily
  when *include_endpoint* is *.true.*.
- *contour_anchor*: manually specify the left endpoint of
  the decay-rate contour. If less than or equal to zero, the code decides the
  left endpoint. This point is, by default, always zero unless a solution using
  the finite-temperature HFB is requested. In this case, the left endpoint is
  by default *E_gs / 2* to avoid both the finite-temperature pole at E = 0 and
  cutting off any strength which lies below our estimate of the ground-state
  energy.
- *poly_order*: the order of the polynomial that is fitted to the phase-space
  integral and then extended to the complex plane (because the phase-space
  integral itself is not complex analytic).
- *poly_fit_si*: the required accuracy of the polynomial fit.
- *poly_fit_points*: the number of points used in the polynomial fit.

### Namelist *SOLVER*

Fine-tuning parameters for the non-linear pnFAM solver.

- *max_iter*: maximum number of iterations for solving a set of FAM equations
- *broyden_history_size*: number of previous iterations considered in Broyden
  mixing
- *convergence_epsilon*: the convergence criterion for the iteration

### Namelist *PARALLEL*

Only read by *pnfam_mpi.x* in the single-file mode, ignored by *pnfam_nompi.x*.

- *nr_adlb_servers*: the number of ADLB servers.  If set to zero, ADLB is
  disabled for the run and points are assigned in round-robin fashion instead
  of dynamic load balancing.  The FINDMAX mode cannot take advantage of the
  dynamic load balancing, so for that mode this parameter is neglected.


Running in a multi-file mode (pnfam_mpi.x only)
-----------------------------------------------

The MPI version of pnFAM can also be run in a multi-file mode, in which several
computations share the same ADLB server, to conserve CPU resources in an MPI environment.
In the multi-file mode, an equally-sized group of processors handles each input file,
although the load balancing only happens within each processor group.

To run *pnfam_mpi.x* in the multi-file mode, you need an additional parameter file
specifying the number of ADLB processes (usually one should be enough, considering the low
communication overhead) and the multiple input files as used in the single-file mode.
This parameter file (see *parallel.in* for an example) consists of a namelist with the
entries

- *nr_adlb_servers*: the number of ADLB servers.
- *nr_parameter_files*: the number of parameter files.

After the namelist, the parameter files should be listed, one file name per line.
Note: Every file in the multi-file mode must run the same type of computation.

To run in the multi-file mode, do not specify any command-line parameters, but direct
the *parallel.in* file to the standard input, e.g.

    mpirun -np 12 ./pnfam_mpi.x < parallel.in


Note on experimental features
-----------------------------

While this release has been stripped of some of the less refined extensions used in the
UNC Chapel Hill research, an experimental finite-temperature mode has been left in the
code. THIS FEATURE IS NOT THOROUGHLY TESTED AND HENCE WE LEAVE IT UNDOCUMENTED. If you
wish to use it, be aware that you may need to debug it and develop it further.
