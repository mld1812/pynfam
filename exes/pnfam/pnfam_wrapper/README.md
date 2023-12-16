The pnFAM Quick Reference
=========================

There are currently four executables in the project: *betadecay.x*, *hfbinfo.x*, and
two pnFAM executables.
- The pnFAM executable named *pnfam_nompi.x* can be compiled to run in serial or with
  OpenMP.
- The one named *pnfam_mpi.x* can be compiled as purely MPI or as a hybrid MPI+OpenMP
  code, and can use ADLB for load-balancing. It parallelizes the number of FAM calculations,
  and has a multi-file mode to process several input files at a time.
- The beta decay executable reads the contour mode binary output files to compute the
  allowed and first forbidden contributions to the rate.
- The hfbinfo executable lists in plain-text the details of any HFBTHO binary output file.
 
Installation
------------

### Dependencies

1. HFBTHO: For input, the pnFAM codes require a binary output file from HFBTHO:
   - HFBTHO v3.00d can be compiled with the "USE_QRPA" flag set to 1 in the Makefile
     to generate the necessary output file. This version is not yet on the CPC
     archive (as of July 09, 2018).
   - HFBTHO v2.00d can be downloaded from the CPC archive but must be modified. The
     required modifications are included with the pnFAM as a patch file (instructions
     below).
2. GSL: Both the MPI and non-MPI versions need to be linked to the GNU Scientific
   Library (GSL).
3. ADLB: The MPI version optionally requires the ADLB load-balancing library,
   available for free from https://www.cs.mtsu.edu/~rbutler/adlb/.
4. Compilers/OS: The programs have been tested to successfully compile with
   sufficiently new gfortran and Intel Fortran compilers on various platforms ranging
   from an Apple MacBook Pro running OS X Yosemite to the Stampede supercomputer in
   Texas Advanced Computing Center (TACC) and Hopper in National Energy Research
   Scientific Computing Center (NERSC).

### Step-by-step installation instructions

1. Download HFBTHO:
   - Version 2.00d is available in the Computer Physics Communications Program
     Library (http://cpc.cs.qub.ac.uk). Copy the file pnFAM file hfbtho.patch to the
     HFBTHO directory and then apply our patch that extends HFBTHO output using the
     command `patch < hfbtho.patch` in the HFBTHO directory.
   - Version 3.00d is located on GitLab with appropriate access (contact
     schunck1@llnl.gov). No modification is required.
3. Modify the HFBTHO Makefile as needed for your platform and compile HFBTHO.
   - Set "USE_QRPA" flag to 1 for HFBTHO v3.00d.
4. (Optional for compiling the MPI version:) Download and compile ADLB
   (https://www.cs.mtsu.edu/~rbutler/adlb/) in a separate directory, following their
   instructions. Then either copy the files _adlbf.h_, _libadlb.a_, and _libfadlb.a_ over
   to the pnFAM directory, or modify the pnFAM makefile so that they are found.
   Note that you must compile pnFAM and ADLB with the same compiler and MPI implementation.
5. Modify the pnFAM makefile for the compiler installed on your platform and so the
   dependency libraries are found. Adjust the options to compile with/without MPI/ADLB/OpenMP.
   Compile the 4 executables via "make all".
6. (Recommended:) Run the Bash script `test_nompi.sh` in the tests/ subdirectory to
   confirm the non-MPI codes were properly compiled. You will likely have to set the
   paths in the beginning of the script to reflect your file organization for the tests to pass.
   Make sure your environment $PATH variable also includes the necessary paths to libraries (e.g. GSL).


Running in a single-file mode
-----------------------------

In a single-file mode, the parameters are read by default from file *default.in*.
A single command line input is accepted for files named differently. The MPI version of
pnFAM parallelizes the number of FAM points to calculate. The input parameters are
interpreted as follows:

### Namelist *GENERAL*

- *fam_mode*: currently the options are "STR" for computing strength functions,
"CONTOUR" for computing closed contours for beta-decay rates, and "FINDMAX" for
locating a (local) maximum of the strength function (useful for finding
resonance peak energies).  No default, must be given explicitly.
- *hfb_input_filename*: the file containing the special binary output from HFBTHO
designed for the pnFAM.
- *fam_output_filename*: the final output file for pnFAM.
- *log_filename*: a file where verbose information about the run should be
printed. If empty, writes to standard output instead.

### Namelist *EXT_FIELD*

- *operator_name*: selected operator to compute.  Currently implemented options
  are:
    + 'F+',   'F-'  : Fermi operator (identity) for beta+ and beta- (J=0)
    + 'GT+',  'GT-' : Gamow-Teller operator (pauli matrix sigma) for beta+ and beta- (J=1)
    + 'RS0+', 'RS0-': position vector (r) and sigma coupled to J=0
    + 'RS1+', 'RS1-': r and sigma coupled to J=1
    + 'RS2+', 'RS2-': r and sigma coupled to J=2
    + 'PS0+', 'PS0-': momentum over i (p/i) dotted with sigma
                      (note: this is not equal to p/i and sigma coupled to J=0)
    + 'R+',   'R-'  : r for beta+ and beta-
    + 'P+',   'P-'  : p/i for beta+ and beta-
 
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
- *half_width*: the imaginary part added to the energies omega. This introduces
  the smearing of the strength function (by a lorentzian with this FWHM).  Note
  that smaller half-width requires more closely spaced points to capture the peaks
  (i.e. smaller *energy_step*) and that in general, the solution of the
  non-linear pnFAM equations takes more iterations when closer to the real axis
  (where the poles of the strength function lie).

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
- *half_width*: the half-width for the strength function, as in STR_PARAMETERS.
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
- *contour_anchor*: specify the left endpoint of the decay-rate contour.
  A sensible default is 0 MeV so that the contour only includes
  positive-energy poles.
- *use_gauleg_contour*: Choose angles *t* as the roots of legendre polynomial so the
  final contour integration for the rate can be done with gauss-legendre quadrature.
  If false, an evenly spaced grid for *t* is used.
- *ratint_psi_gl_points*: A preliminary rate for the operator is computed by the pnfam
  executable. The phase space integrals needed for this are computed by default with
  gauss-legendre quadrature. This parameter is the number of quadrature points to use.
- *ratint_points*: The phase space integral as a function of the upper integration limit
  is calculated at points along the real axis. By default these points are roots of a
  chebychev polynomial, which are then interpolated with a rational-function so it can
  be analytically continued. This parameter is the number of phase space integrals to
  compute for the interpolation. Note that the beta-decay executable has more options
  for computing the phase space, this is just a default method for a preliminary calc.

### Namelist *SOLVER*

Fine-tuning parameters for the non-linear pnFAM solver.

- *max_iter*: maximum number of iterations for solving a set of FAM equations
- *broyden_history_size*: number of previous iterations considered in Broyden
  mixing
- *convergence_epsilon*: the convergence criterion for the iteration


Running in a multi-file mode (pnfam_mpi.x only)
-----------------------------------------------

The MPI version of pnFAM can also be run in a multi-file mode. If ADLB is used,
several computations share the same ADLB server, to conserve CPU resources in an MPI
environment. In the multi-file mode, an equally-sized group of processors handles
each input file, although the load balancing only happens within each processor
group.

To run *pnfam_mpi.x* in the multi-file mode, you need an additional parameter file
specifying the number of ADLB processes (usually one should be enough, considering
the low communication overhead) and the multiple input files as used in the single-
file mode. This parameter file (see *parallel.in* for an example) consists of a
namelist with the entries:

- *nr_adlb_servers*: the number of ADLB servers. If set to zero, ADLB is disabled for
  the run and points are assigned in round-robin fashion instead of dynamic load
  balancing. The FINDMAX mode cannot take advantage of the dynamic load balancing, so
  for that mode this parameter is neglected.
- *nr_parameter_files*: the number of parameter input files.

After the namelist, the parameter files should be listed, one file name per line.
Note: Every file in the multi-file mode must run the same type of computation.

To run in the multi-file mode, do not specify any command-line parameters, but direct
the *parallel.in* file to the standard input, e.g.

    mpirun -np 12 ./pnfam_mpi.x < parallel.in


Combining contour integrals to a total beta-decay half-life (betadecay.x)
-------------------------------------------------------------------------

While the total allowed half-lives can be easily obtained by manually adding up the
Fermi and Gamow-Teller contributions, the case where first-forbidden channels are
included is more complicated. To sum up the different contributions correctly, the
auxiliary program *betadecay.x* should be used. This program also allows the use of
different numerical methods for computing the phase space and contour integrals.

The input parameters for *betadecay.x* are placed in the file *beta.in*, and they are
interpreted as follows:

### Namelist *BETA_PARAMETERS*

- *hfbtho_solution_file*: the HFB solution file produced by HFBTHO.

If left empty, these reamining parameters are read from the HFB solution file. In
other words, specifying these in *beta.in* can be used to override the values in the
HFB solution file.

- *z_init*: the charge number of the mother nucleus.
- *a_init*: the mass number of the mother and the daughter nucleus.
- *q_value*: the Q value of the decay.
- *gs_energy*: the energy of the lowest QRPA excitation.
- *ga_effective*: the effective axial-vector coupling constant (default: 1.0,
   representing the quenching of the axial-vector current in nuclei).
- *output_file*: the file name for a more machine-readable output file; if omitted,
   the output is stored in "beta.out".

### Namelist *PHASE_SPACE*

- *quadrature_type*: the type of numerical quadrature used to compute the contour
  integration. Options are "GAUSS","SIMPSON", or "TRAP" for gauss-legendre,
  simpson's 3/8, and trapezoidal rule, respectively. Note that GAUSS can only be used
  if pnfam was run with use_gauleg_contour = true. No default, must be specified explicitly.
- *use_polyfit*: if true, the phase space integrals are analytically continued using
   a least squares polynomial fit (with GSL) to points along the real axis. If
   false, the same is accomplished by interpolation via a rational function.
- *polynomial_order*: If use_polyfit = true, the order of the polynomial fitted to the
   phase-space integrals.
- *polynomial_fit_points*: If use polyfit = true, the number of points used in the
   polynomial fit.
- *polynomial_fit_si*: If use_polyfit = true, the required accuracy of each phase
   space integral, calculated by simpson's 3/8 with mesh refinement until si achieved.
   This method is used only to maintain consistency with previous versions of pnFAM.
- *ratint_psi_gl_points*: If use_polyfit = false, each phase space integral is
   calculated with gauss-legendre quadrature with this many points.
- *ratint_points*: If use_polyfit = false, number of phase space integrals to use for
   rational function interpolation. These are computed at chebychev roots along the
   real axis.

### Namelist *INPUT_FILES*

If any of the following file names are empty, the corresponding (K, Parity)
contribution to the decay rate is simply omitted.  E.g. omitting *k1_rs1* will result
in the exclusion of the K=1 first forbidden decay rate.

- *k0_f*: the pnFAM binary solution file containing the pnFAM solutions at the
   contour points for the K=0 component of the Fermi operator.
- *k0_rs0*: the same for the RS0- operator.
- *k0_ps0*: the same for the PS0- operator.
- *k0_gt* and *k1_gt*: the same for the K=0 and K=1 components of the Gamow-Teller
   operator, respectively.
- *k0_p* and *k1_p*: the same for the K=0 and K=1 components of the P- operator,
   respectively.
- *k0_r* and *k1_r*: the same for the K=0 and K=1 components of the R- operator,
   respectively.
- *k0_rs1* and *k1_rs1*: the same for the K=0 and K=1 components of the RS1-
   operator, respectively.
- *k0_rs2*, *k1_rs2*, and *k2_rs2*: the same for the K=0, K=1, and K=2 components of
   the RS2- operator, respectively.


Note on experimental features
-----------------------------

While this release has been stripped of some of the less refined extensions used in
the UNC Chapel Hill research, an experimental finite-temperature mode has been left
in the code. THIS FEATURE IS NOT THOROUGHLY TESTED AND HENCE WE LEAVE IT
UNDOCUMENTED. The finite-temperature mode is only activated when the HFB ground state
was obtained from HFBTHO with T not equal to zero. If you wish to use this extension,
be aware that you may need to debug it and develop it further.
