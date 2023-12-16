The pnFAM Quick Reference - Contour Executable
=========================

pnFAM now has the ability to include two-body beta-decay operators 
coming from the two-body nuclear current. These are treated in mean-field
theory and act as corrections to the typical one-body operators coming
from the one-body current. As such, there are no "new" operators for
the two-body current. Instead, there is a two-body current mode, which
specifies whether to compute to one-body, two-body, or one-plus-two-body
contribution for the standard one-body operators described in the main
pnfam documentation.

Status of Implementation
------------------------

As of June 2021, only the two-body mean-field contribution to the
Gamow-Teller operator is implmeneted.

Two-body Gamow-Teller pairing-field is in progress. pnFAM matrices have
only one block per row and column, but the full two-body external field
with pairing would break this structure in the QP basis. To implement
pairing, we must extend the block matrix type to handle completely general
block structures. This can be accomplished by replacing the r2c and c2r
methods with lookups in two indpendent arrays, one for block rows and one
for block columns.

Forbidden operators have not been considered.

Implementation Overview
-----------------------

A two-body external field appears in linear response theory in the exact same
way as the two-body hamiltonian. Therefore including it requires computing
two-body matrix elements of the external field and then contracting those
matrix elements with the one-body density matrices to obtain a one-body operator.

### The Problem:
The spatial pieces of the two-body matrix elements of the current are tricky because
they contain Yukawas and derivatives of Yukawas. Yukawas are problematic because
for two-body matrix elements in:

1. coordinate space, the finite range requires a very large integration mesh
2. configuration space, they are not separable in coordinates and therefore require
   and unreasonable amount of nested loops.

### Using Gaussians
To combat this, we approximate the Yukawas as a sum of gaussians. Yukawas are singular at r=0
and gaussians are not, their matrix elements contain an extra factor of r^2, and so the matrix
element *integrand* of Gaussians and Yukawas can agree quite well. CPC 180 (2009) 1013-1040
desribes how two-body matrix elements of Gaussians can be computed in the axial HO basis using
analytic functions of the HO quantum numbers. Derivatives of Yukawas are treated by integrating
by parts to act the derivatives on the HO wavefunctions. These in turn can be expressed as sums
of other HO wavefunctions, allowing the use of the analytic matrix element expressions everywhere.

### Cartesian Expansion
The z-component of Gaussian matrix elements is simplest, because it involves 1D cartesian HO
wavefunctions. The (radial,angular)-component is more complicated because the two coordinates are
coupled. As described in the CPC paper, one can expand the (radial,angular) part of the HO wavefunctions
in 1D cartesian HO wavefucntions. A Gaussian is also trivially separable in all 3 cartesian coordinates.
Using this, we can break the matrix element V(r,phi,z) into V(r,phi) = C(x,y)*V(x)*V(y) and V(z).
This allows us to only need 1D cartesian matrix elements of a Gaussian and derivatives of 1D cartesian
HO wavefunctions.

### Computational Scaling
The best method for computing the two body contributions is to work in configuration space and take
advantage of separability in the interaction. The following discuss this scalability of this approach
in terms of the number of oscillator shells N.

Yukawas are not separable, so to contract their matrix elements with the 1-body densities we need
6 loops over the 3 quantum numbers to populate the 2 indices of the mean field, and 6 loops to contract
the remaining 2 indices of the two-body matrix elements with the density, which is O(N^12).

Gaussians and axial HO wavefunctions can both separate into f(r,phi)*g(z). The contraction can then
separate into 2 sets of nested loops, one treating just the (r,phi) quantum numbers and the other
treating just the (z) quantum numbers. We still need 6 loops over all three quantum numbers to populate
the 2 indices of the mean field, but the contraction ove the remaining 2 indices is split into a
contraction over just (r,phi) quantum numbers, which is then O(N^10), followed by a contraction over
just (z) quantum numbers, which is O(N^8).

For comparison, computing matrix elements for the 1-body external field contains 2 loops over the whole
basis. The basis size is N*(N+1)(N+2)/6 ~ N^3, so the 1-body extfield is O(N^6).

### Solution summary:

1. We express the current in spherical components so it aligns with the existing
   implementation. For example, GT operator is split into K=+1,0,-1, each of these components
   will have a 2BC correction.
2. We fit a small set of gaussians to the Yukawa and use integration by parts for derivatives
   to express the current matrix elements entirely in terms of matrix elements of gaussians.
   We then borrow existing gaussian matrix element routines from HFBTHO.
3. We contract the two-body matrix elements with the one-body density matrices to construct the
   two-body external field, using separability of the interaction to make the computation feasible.

Code structure
-----------------------
Several modules were added to explicitly handle the two-body current:

- pnfam_spatial_mtxels.f90
- pnfam_type_extfield_2bc.f90
- pnfam_extfield_2bc.f90

pnfam_spatial_mtxels.f90 contains the methods for computing two body matrix elements
of gaussians and derivatives of gaussians. pnfam_type_extfield_2bc.f90 contains methods
to compute and store two-body matrix elements of the current operator. pnfam_extfield_2bc.f90
computes the contraction of the current with the one body density matrices.

The external field setup is moved from pnfam_setup.f90 to pnfam_solver.f90, where the two-body
current can be added to the usual one-body current external field, or computed by itself, depending
on the mode. pnFAM is then run as usual using the resulting external field matrix.

Storage
------------------------
The contraction of the two-body current with the density matrices is extremely expensive.
Presently, highly parallelized calculations using OpenMP for the for-loops takes about
1.5-2.5 hours for a 16 shell calculation (per K-component). Since the resources needed for
the external field computation are no longer negligible, the contracted two-body external field
is stored in a file with extension .tbc. Future calculations for the same nucleus can skip the
computation by reading from this file.

The contact part of the two-body current is exactly proportional to the density times the usual
gamow-teller operator, and can be computed cheaply using the one-body external field routines.
Therefore, to provide the most flexibility, only the finite-range part of the current is calculated
in pnfam_extfield_2bc.f90 and stored in the .tbc file (c3, c4 and P terms). The contact term is 
re-computed using the one-body routine every execution.

The finite range mean field matrix is split into 6 pieces: a c3, c4, and P piece for direct
and exchange. They are stored with the LECs stripped off, such that when the .tbc file is read
any new set of LECs from the namelist can be applied.

Approximation schemes
-----------------------
We have also included some approximations to the two body currents to compare with the full
calcualtion. The first scheme uses nuclear matter to reduce the current to an effective one-body
operator. We have implemented the result for symmetric and isospin-asymmetric nuclear matter, both
of which result in density dependent functions mulitplying the usual gamow-teller operator. The
expressions are implmenented for zero momentum transfer (valid for beta decay) and zero total
momentum. More general expressions exist for non-zero total momentum but are not implemented.

We have also added the density matrix expansion (DME) for the current, which expresses
the two-body current as an expansion in derivatives of the one-body density matrices (and so is
itself an effective one body operator). So far we have implemented the following terms:

- Direct: LO(only contact is non-zero), NLO(is zero), N2LO, N3LO(is zero), N4LO
- Exchange: LO

Parallelization
-----------------------
The openMP computation of the full FAM 2BC requires openMP to run in a reasonable time scale.
For a slurm scheduler, pnFAM should be run with 1 MPI task and many openMP threads for that task.
Avoiding multi-threading should provide optimal performance. As the number of oscillator shells
becomes large, the openMP stack size must be increased to avoid segmentation faults. An example
of slurm settings for a 16 shell calculation is given below:
"""
#!/bin/bash
#SBATCH -p example_queue
#SBATCH -J example_jobname
#SBATCH -o slurm_%A.out
#SBATCH -t 0-04:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=44
#SBATCH --threads-per-core=1
#SBATCH --hint=nomultithread
export OMP_NUM_THREADS=44
export OMP_STACKSIZE=256m
export OMP_DYNAMIC=false
export OMP_PROC_BIND=true
export OMP_WAIT_POLICY=active
ulimit -s unlimited
./pnfam_main.x GT-K0.in
"""
