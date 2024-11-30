The pnFAM Quick Reference
=========================

The program pnFAM solves the Finite Amplitude Method equations for a single point
on the QRPA strength function for one of a set of charge-chaning external field operators.
It utilizes interactions derived from generalized energy density functionals and is
self-consistent with Hartree-Fock-Bogoliubov ground states computed with the same functional.
The HFB solver HFBTHO is integrated into pnFAM, which starts the computation from the HFBTHO
output binary file. The following symmetries are assumed in the pnFAM solver:
- Isospin symmetry in the HFB ground state (no proton-neutron mixing)
- Time-reversal symmetry in the HFB ground state (equal filling for odd states)
- Axial symmetry (Matrices come in blocks with well defined intrinsic Z-projection of Ang. Mom.)
pnFAM assumes the Bogoliubov transformation provided from the HFB solution makes use of the
Bloch-Messiah theorem, such that the matrices in pnFAM take on a particular block structure:
for every block-row, there is only one non-zero block-column.

Installation
------------

### Dependencies

1. HFBTHO: For input, the pnFAM codes require a binary output file from HFBTHO:
   - HFBTHO v4.00 (included in this repository) can be compiled to produce libhfbtho.a
     static library. This is linked when compiling pnFAM v2.00, allowing it to use
     "hfbtho_output.wel" as the pnFAM input file. The "hfbtho_NAMELIST.dat" is also required.
   - HFBTHO v3.00 can be downloaded from the CPC archive and compiled with the "USE_QRPA"
     flag set to 1 in the Makefile. This generates the "solution.hfb" file which is
     the input file used by pnFAM v1.00.
2. Compilers/OS: The programs have been tested to successfully compile with
   sufficiently new gfortran and Intel Fortran compilers on various platforms ranging
   from an Apple MacBook Pro running OS X Yosemite to the Stampede supercomputer in
   Texas Advanced Computing Center (TACC) and Hopper in National Energy Research
   Scientific Computing Center (NERSC).

### Step-by-step installation instructions

1. Download HFBTHO:
   - Version 4.00 is provided in this repository
   - Version 3.00 is available in the Computer Physics Communications Program
     Library (http://cpc.cs.qub.ac.uk).
2. Modify the HFBTHO Makefile as needed for your platform and compile HFBTHO.
   - Compile the library libhfbtho.a for HFBTHO v4.00.
   - Set "USE_QRPA" flag to 1 for HFBTHO v3.00.
3. Modify the pnFAM makefile for the compiler installed on your platform and so the
   dependency libraries are found.
   - For pnfam v2.00, the compilation settings should match those used to compile
     libhfbtho.a (e.g. same compiler, both with openmp disabled, etc.)

Running the pnFAM program
-----------------------------

The pnFAM input parameters are read by default from file *pnfam_NAMELIST.dat*.
A single command line argument is accepted for files named differently.

### Namelist *GENERAL*

- *fam_output_filename*: String, to be used as the base name for output filenames.
  The text output file will have extension '.dat' and the binary will have '.bin'.
  If blank or empty string, no output files are written, and the following values are
  imposed for the inputs below: print_stdout=T, use_fam_storage=0
- *print_stdout*: Boolean, indicating whether the program should print output to stdout.
- *use_fam_storage*: Integer specifying use of stored pnfam solution.
  For  0, pnFAM will neither read nor write the solution
  For  1, pnFAM will read the soltuion (if it exists) but will not write the solution
  For -1, pnFAM will read the soltuion (if it exists) and will write the solution
  The fam solution file name is specified via fam_output_filename. For -1 an existing
  file will be overwritten with the new solution.
- *real_eqrpa*: Real, for the real part of the energy at which the strength function
  is computed.
- *imag_eqrpa*: Real, for the imaginary part of the energy at which the strength
  function is computed. This is equal to half the full-width-at-half-max of a Lorentzian
  function which smears poles of the strength function in the complex plane.

### Namelist *EXT_FIELD*

- *beta_type*: String, to indicate the type of decay/external field operator. "-" for
  beta-minus decay (isospin lowering operator) or "+" for beta-plus decay and
  electron-capture (isospin raising operator).
- *operator_name*: String, for the selected operator multiplying the isospin ladder
  operator defined through *beta_type*. The operators, coming from various beta decay
  multipoles, are implemented by spherical K-component. This means we re-couple
  products of operators where applicable onto definite J, and implement the components
  K=-J,...,+J separately. Implemented options are:
    + 'F'  : Fermi operator (the identity)              (J=0, parity=+)
    + 'GT' : Gamow-Teller operator (pauli matrix sigma) (J=1, parity=+)
    + 'R'  : position vector (r)                        (J=1, parity=-)
    + 'P'  : momentum divided by i (p/i)                (J=1, parity=-)
    + 'PS0': p/i dotted with sigma                      (J=1, parity=-)
    + 'RS0': r and sigma coupled to J=0                 (J=0, parity=-)
    + 'RS1': r and sigma coupled to J=1                 (J=1, parity=-)
    + 'RS2': r and sigma coupled to J=2                 (J=2, parity=-)
      (Note: dot product in PS0 is not equal to tensor product coupled to zero.
       For vector operators, A.B = -sqrt(3) [AB]_00. See "From Nucleons to Nucleus"
       by J. Suhonen Ch. 2 for more details).
- *operator_k*: Integer, indicating the spherical component of the operator defined by
  *operator_name* to compute (K=-J,...,+J). This is the the amount by which the operator
  changes intrinsic Z-projection of angular momentum. Since time-reversal is conserved,
  negative values yield the same result as positive ones.
- *compute_crossterms*: Boolean, indicating whether to compute, in addition to the
  strength, the cross-terms of the selected operator needed to evaluate non-unique
  forbidden decays. These are obtained by tracing the FAM result for operator O1 with
  matrix elements of operator O2, giving result proportional to <O1><O2>.
- *two_body_current_mode*: Integer, 3 digits indicating treatment of two body weak currents.
  Nuclear matter LDA approximations do not include pairing, so the only non-valid combination
  of the digits below is if 2nd digit is not 1, the 3rd digit cannot be 2 or 3.
    + To not use two body currents at all, supply a value of 0, otherwise supply a 6 digit integer.
    + 1st digit: 1=1BC+2BC, 2=2BC only
    + 2nd digit: 1=Full pnFAM 2BC, 2=SNM+LDA 2BC, 3=ASNM+LDA 2BC, 4=DME(exc), 5=DME(exc)+FAM(dir)
    + 3rd digit: 1=No pairing (Gamma only), 2=Pairing only (Delta only), 3=Gamma and Delta
    + 4th, 5th, 6th digits: Determines which terms to use (4th: GT, 5th: P (vector), 6th: PS0 (axial charge)).
    + For the 4th digit: 1 for 2-body GT, 2 for 2-body GT + spin-dipole. 5th and 6th digits: 1 for 2-body fermi gas, 2 for DME.
- *two_body_current_usep*: Boolean, if true uses P-terms in 2BC
- *two_body_current_lecs*: Real, three values defining the low energy couplings c3, c4, cd.
  The LECs supplied should be dimensionless as in PRC 67 055206 (2003) Eq.12. cd = d1 + 2\*d2.

### Namelist *INTERACTION*

- *interaction_name*: String, for the selection of a built-in interaction or an
  external coupling constant file. An external file 'pnfam_CUSTOM_EDF.dat' would be
  selected by specifying 'FILE:pnfam_CUSTOM_EDF.dat' here (the file name can be user defined).
  The residual interaction can be turned off completely by specifying 'NONE'.
- *require_self_consistency*: Boolean, indicating whether to proceed with the calculation
  if the coupling constants of the residual interaction are not consistent with those
  of the HFB solution being used. If .true. the program will abort if the
  couplings are not the same (within tolerance). Used to guard against unintentional
  breaking of self consistency of the interaction.
- *require_gauge_invariance*: Boolean, if .true. the program (1) adjusts the time-odd
  coupling constants CT and CF to obey the gauge relations linking them to the
  J^2 couplings and then (2) aborts if any coupling constants do not obey their
  individual gauge relations. Used to guard against unintentional breaking of
  the gauge invariance.
- *force_j2_terms*: Boolean, if .true. the J^2 terms are included in the computation,
  even if they shouldn't be ordinarily (e.g. forces such as SLy4 or SkM\* where
  the J^2 terms are neglected in the HFB computation).  If .false., the
  J^2 terms are chosen to match what was used in the HFB computation (e.g. for
  SLy4 there are no J^2 terms, but for SLy5 they are included).
- *vpair_t0*: overrides the default isoscalar pairing strength
- *vpair_t1*: overrides the default isovector pairing strength
- *override_cs0*:  overrides the default spin-density coupling (density-independent part)
- *override_csr*:  overrides the default spin-density coupling (density-dependent part)
- *override_cds*:  same for the laplacian s term
- *override_ct*:   same for the T term
- *override_cf*:   same for the F term
- *override_cj*:   same for the j term
- *override_cgs*:  same for the div s term
- *override_csdj*: same for the curl j term
(Note: the default couplings are derived from the Skyrme interaction parameters t,x)

### Namelist *SOLVER*

- *max_iter*: Integer, maximum number of iterations for solving the FAM equations
- *convergence_epsilon*: Real, the convergence criterion for the iteration, defined as the
  maximal change of any matrix element in the FAM amplitudes (P,X,Y,Q) from one iteration
  to the next.
- *broyden_history_size*: Integer, number of previous iterations retained in Broyden mixer.
  Note: this is memory intensive, creating an array with size nbroy, where
  nbroy = (broyden_history_size)\*2\*(number fam amplitudes)\*(non-trival fam matrix elements)
  The number of FAM amplitudes is 2 (X,Y), or 4 if using statistical FAM (P,X,Y,Q)
- *energy_shift_prot*: Real, an amount by which to shift all proton QP energies
- *energy_shift_neut*: Real, an amount by which to shift all neutron QP energies
  Note: shifting the QP energies spoils the QRPA eigenvalue problem and will change the
  the magnitude of the strength, in addition to the location of the poles.
- *quench_residual_int*: Real, an amount by which to multiply the residual interaction
  matrix. A value of 0.0 is the same as using *interaction_name* = 'NONE', and a
  value of 1.0 uses the unaltered residual interaction.

Note on statistical FAM calculations
---------------------------------------

### Equal Filling Approximation

pnFAM has the ability to calculate nuclei with odd particle number in the equal filling
approximation. This is a time-reversal invariant statistical ensemble in which the odd
quasiparticle and it's degenerate partners are the only occupied states, and they are
occupied with equal probability. In the present case of axially symmetric systems, this
is just two states, with occupations f=1/2. These calculations are explored in detail in
PRC 102, 034326 (2020). That reference notes several problems with EFA calculations.

1. Positive strength transitions may appear at negative QRPA energies
2. Negative strength transitions may appear at positive QRPA energies
3. Both 1. and 2. may appear in the same strength function ("overlapping transitions")
4. Imaginary poles may be present, which are difficult to identify without a full 2D contour

Note that while in some cases it might appear that 1. and 2. are occuring on their own,
i.e., the transition at E does not appear to have a partner at -E, this is not the case.
The QRPA transitions are ALWAYS symmetric about E=0. Observing this situation just means
that the strength/residue of the "missing" transition is negligibly small, but can still
be observed by increasing the resolution (significantly).

*I believe that for the stability condition for statistical QRPA to be met, the quasiparticle
occupations must be a part of the minimization. This is because the P and Q terms of the response
account for transitions AMONG occupied states. If the occupations are not part of the minimization,
some of these transitions may be to states that were not considered during the minization and that
have lower energy. The EFA does not inlcude occupations in the minimization explicitly because they
are manually set to 1/2 for the chosen odd state and its time-reversed partner. However, they are
included indirectly by running over all possible states in the blocking procedure and choosing the
configuration that gives the lowest energy. This MIGHT suffice for the like-particle QRPA, as long
as EFA ground state has the lowest energy QPs blocked (which might NOT be the case for odd-odd
nuclei due to Coulomb), but in the charge changing channel the blocking procedure NEVER explores
all possible n-QP configurations in the EFA ensemble that can be included by the P and Q terms.
Furthermore, while in principle the problems 1-3 listed above are stable transitions to states of
lower energy (because the energy is still real) I believe they are still problematic. They are
simply "unstable" transitions that are SO unstable that the QRPA residual interaction is not large
enough to make it imaginary. The above claim remains to be rigorously proven, but use the EFA with care.*

### Finite Temperatures

pnFAM has the ability to calculate strength functions at finite temperatures. This feature
is automatically turned on if the HFBTHO ground state is computed at finite temperature.
pnFAM computes the "bare" finite temperature response, which includes a sum over not just final states,
but also initial states. Each term in the sum is weighted with an ensemble weight that depends on the energy
of the initial state, Pi=exp(-Ei/kT). This means the QRPA response contains both transitions from E1 to E2
as well as E2 to E1, but the come with different weights.

A problem with the bare response function though is that the E1 to E2 transition provides a positive contribution
to the strength at E=+(E1-E2), while the E2 to E1 transition provides a negative contribution at the SAME energy
E=-(E2-E1), coming from the reverse process. The latter transition however will be weighted with a different value
P2=exp(-E2/kT) than the former’s weight P1=exp(-E1/kT). The weights cause the net strength at E>0 to always be
positive, but reduced. A similar thing happens at negative energy, where the net strength is always negative.

Therefore, to get the finite temperature strength from the bare response we must remove contributions from the
reverse process, because they overlap with the transitions from the forward process. This is done via detailed balance,
and amounts to multiplying the bare response by a prefactor 1/(1-exp(-E/kT)), where E is the QRPA energy. This
prefactor increases that reduced strength at E>0 by the necessary amount, and makes the strength at E<0 positive by the
correct amount. We should apply this prefactor, then use strength at both positive and negative energies. The
strength at E<0 corresponds to excited parents de-exciting, and has a larger Q value than the strength at the same |E|>0,
which is the typical transition where the parent ground state is excited. The prefactor dies to 0 at negative energies,
and approaches 1 at positive energies, so it has a temperature-dependent range of importance in the region near E=0.

This prefactor is NOT included in the pnFAM result. pnFAM only calculates the bare response. The reason is that
the energy inside the prefactor depends on the type of energy contour, while pnFAM computes only a single energy
point and has no concept of the contour. For example, we pretend strength calculations are on the real axis,
and so we should use a real energy in the prefactor. For contour calculations, however, we should use the full
complex energy in the prefactor. Note that the prefactor contains a pole at E=0, which must be treated with care.
Also, whether detailed balance can be applied to beta-decay is not strictly justified. We can assume detailed
balance under the approximation where there is zero energy loss due to the neutrino.
