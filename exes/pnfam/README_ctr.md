The pnFAM Quick Reference - Contour Executable
=========================

The pnFAM repository contains, in addition to pnfam_main.x, and addtional
executable named contour_main.x. This program wraps the main pnfam program, and
handles high-level functionality for computing QRPA strength functions along
discrete energy contours for one or more beta decay operators.

Installation
------------

Installation and dependencies are identical to those listed in the README.md
file for pnfam_main.x

1. (Recommended:) Run the Bash script `test_nompi.sh` in the tests/pnfam2_serial
   subdirectory to confirm the program was properly compiled. You will likely have to
   set the paths in the beginning of the test script to reflect your file organization.

Running the contour program
-----------------------------

The contour input parameters are read by default from file *pnfam_CONTOUR.dat*.
A single command line argument is accepted for files named differently.

### Namelist *CTR_GENERAL*

- *fam_mode*: currently the options are "STR", "CONTOUR", and "FINDMAX".
  "STR" computes along a contour with fixed half-width (for strength functions).
  "CONTOUR" computes along a closed circular contour (for summed strength).
  "FINDMAX" is similar to "STR" mode, but is used for locating local maxima in the
  strength function.
- *fam_input_filename*: the namelist used by pnfam_main.x. The majority of the settings
  for the run are defined in this file, with a few (such as the energy and operator)
  being overridden by the contour functionality. If fam_output_filename is not empty,
  it is overridden with "Operator//beta_type//K//operator_k". If it is empty, it behaves
  the same as with pnfam_main for each fam point, however, a summary file with the base
  file name above is always written.

### Namelist *CTR_EXTFIELD*

- *operator_group*: J-parity label for a group of operators (consult README.md for
  list of operators, Js and parities). This input is informational only and is not
  meant to be changed.
- *operator_active*: List defining which operator groups will be computed, 0 will not
  be computed, 1 will be computed. If every member is 0, only the operator defined in
  the file *fam_input_filename* will be computed.
  NOTE: The operators will be computed with sign(K)=sign(operator_k) as defined in
  the file specified by fam_input_filename.

### Namelist *STR_PARAMETERS*

This namelist is only read if *fam_mode* is set to "STR".
It contains the settings for defining the energy contour for this mode.

- *energy_start*: the starting value for the real part of the QRPA energy
- *energy_step*: the spacing of (real) energy values
- *nr_points*: the number of energy values to compute.
- *half_width*: the fixed imaginary party of the contour energies
  This introduces a smearing of the strength function by a lorentzian of FWHM/2 equal to this.
  Note that smaller half-width requires more closely spaced points to capture the peaks
  (i.e. smaller *energy_step*) and that in general, the solution of the
  non-linear pnFAM equations takes more iterations when closer to the real axis
  (where the poles of the strength function lie).


### Namelist *CONTOUR_PARAMETERS - NOT IMPLEMENTED*

This namelist is only read if *fam_mode* is set to "CONTOUR".
It contains the settings for defining the energy contour for this mode.


### Namelist *FINDMAX_PARAMETERS - NOT IMPLEMENTED*

This namelist is only read if *fam_mode* is set to "FINDMAX".
In the FINDMAX mode, the code operates similarly to STR mode. However, it computes
strength by starting at some intial (real) EQRPA, and increasing this value until
it finds an interval where the sign of the derivative of the resulting strength function
flips. It then bisects until the interval where the sign flip occurs is as small as requested.

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


### Namelist *BETADECAY - NOT IMPLEMENTED*

This namelist input is only applicable to STR and CONTOUR modes.
If active, once the strength calculations are complete, the code computes
a complex-analytic approximation to the beta decay phase space integrals,
constructs a phase-space-weighted shape factor from the computed strength
functions, and integrates this to obtain a beta decay rate.


Running in parallel mode - NOT IMPLEMENTED
-----------------------------------------------

The contour program can be compiled with the use_mpi flag set to 1 to allow
parallel computation. The code generates a list of tasks of size nr_points*nr_operators
(nr_operators is at most 14) which are distributed round-robin to mpi processes
(ADLB can be used?). Once each point is computed, the strength (and cross terms) are
stored in the task. After all tasks are completed, summary output files are written
by rank 0.

For example, to execute a parallel calculation with 12 mpi processes, run:

    mpirun -np 12 ./contour_main.x
