#!/usr/bin/env python
"""
Run a short pnFAM calculation (PS0, K=0) for a broad spread of odd-A nuclei,
to produce the blocked-nucleon spin-density dumps (_spin_neu.dat / _spin_pro.dat)
used by spin_term_size_test.py / spin_term_survey.py / spin_current_q.py.

Modeled on run_pynfam5.py: every hfb/fam/ctr/psi setting is copied from that
script unchanged except proton_number, neutron_number, and which species is
blocked, which vary per nucleus in NUCLEI below.

Each nucleus gets its own output directory (e.g. '163Gd_spindensity_test'), grouped
under the single SCRATCH_DIR subfolder below rather than scattered directly in
tests/, matching the existing 189Gd/163Gd naming convention that spin_term_survey.py's
nucleus-label regex already parses. This means calling pynfam_mpi_calc once per
nucleus in a plain loop (identically on every MPI rank) rather than using pynfam's
native tuple-batch mode, which would put every nucleus under one shared directory
with non-descriptive 000000/000001/... labels instead.

After each nucleus finishes, the fam_soln/fam_meta/<OP-K>.tar produced by pynfam
is extracted in place (rank 0 only), so the raw *_spin_*.dat files described above
end up on disk exactly where the analysis scripts expect them:
    <SCRATCH_DIR>/<outputs>/000000/fam_soln/fam_meta/PS0-K0/000000/PS0-K0_spin_{neu,pro}.dat

Edit NUCLEI to add/remove nuclei. Each entry is (Z, N); the odd one of the two is
blocked automatically (both must not be odd -- odd-odd nuclei are skipped, since
the point is a single, unambiguous blocked-nucleon spin density).
"""
from __future__ import unicode_literals
import os
import glob
import tarfile
from copy import deepcopy

from pynfam import pynfam_mpi_calc

# ------------------------------------------------------------------------------
# All nuclei are grouped under one subfolder of tests/, rather than dumping
# 8 more top-level directories directly into tests/.
# ------------------------------------------------------------------------------
SCRATCH_DIR = './tests/spinsurvey'

# ------------------------------------------------------------------------------
# MPI rank (0 in serial / non-MPI launches) -- only rank 0 extracts tar files.
# ------------------------------------------------------------------------------
try:
    from mpi4py import MPI
    _RANK = MPI.COMM_WORLD.Get_rank()
except ImportError:
    _RANK = 0

# ------------------------------------------------------------------------------
# Element symbols, Z=1..100 (for naming output directories / printouts).
# ------------------------------------------------------------------------------
ELEMENTS = [
    None, 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al',
    'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co',
    'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I',
    'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
    'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au',
    'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U',
    'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
]

# ------------------------------------------------------------------------------
# Survey list: (Z, N). Spans light<->heavy, near-shell<->mid-shell/deformed,
# and odd-N<->odd-Z, to check whether the spin-density / spin-current DME size
# (see spin_term_size_test.py, spin_current_q.py) is generic or case-specific.
# Edit freely -- this is a physics choice, not a code constraint.
# ------------------------------------------------------------------------------
NUCLEI = [
    (8,  9),    # 17O   light,  odd-N, one neutron above the N=8  shell closure
    (9,  10),   # 19F   light,  odd-Z, one proton  above the Z=8  shell closure
    (50, 83),   # 133Sn medium, odd-N, one neutron above the N=82 shell closure
    (64, 99),   # 163Gd heavy,  odd-N, mid-shell / deformed        (already run)
    (64, 125),  # 189Gd heavy,  odd-N, near the N=126 shell closure (already run)
    (65, 94),   # 159Tb heavy,  odd-Z, mid-shell / deformed
    (83, 126),  # 209Bi v.heavy,odd-Z, one proton above doubly-magic 208Pb core
    (94, 145),  # 239Pu v.heavy,odd-N, deformed actinide
]


def make_inputs(Z, N, outputs_name):
    """pynfam_inputs / override_settings for one nucleus, copied from
    run_pynfam5.py verbatim except outputs name, Z, N, and blocking species."""

    if (N % 2 == 0) == (Z % 2 == 0):
        raise ValueError(f"Z={Z}, N={N} is not odd-A (need exactly one odd); "
                          f"odd-odd/even-even are not supported by this survey.")
    n_block = [1, 0, 0, 0, 0] if N % 2 else [0, 0, 0, 0, 0]
    p_block = [1, 0, 0, 0, 0] if Z % 2 else [0, 0, 0, 0, 0]

    pynfam_inputs = {
        'directories': {
            'outputs': outputs_name,
            'exes': './exes',
            'scratch': SCRATCH_DIR
        },
        'nr_parallel_calcs': None,
        'rerun_mode': 0,
        'hfb_mode': {
            'gs_def_scan': (-2, (-0.2, 0.0, 0.2)),
            'dripline_mode': 0,
            'ignore_nonconv': 1
        },
        'fam_mode': {
            'fam_contour': 'CIRCLE',
            'beta_type': '-',
            'fam_ops': ('PS0', 0)
        }
    }

    override_settings = {
        'hfb': {
            'proton_number': Z,
            'neutron_number': N,
            'number_of_shells': 16,
            'functional': 'SKOP',
            'number_iterations': 600,
            'accuracy': 1e-08,
            'proton_blocking': p_block,
            'neutron_blocking': n_block,
        },
        'ctr': {
            'nr_points': 2
        },
        'fam': {
            'use_fam_storage': 0,
            'interaction_name': 'SKOP',
            'force_j2_terms': False,
            'convergence_epsilon': 1e-07,
            'require_self_consistency': True,
            'require_gauge_invariance': True,
            'vpair_t1': None,
            'max_iter': 300,
            'broyden_history_size': 200,
        },
        'psi': {
            'GA': -1.0,
            'UseNewShapeFactors': False
        }
    }

    return pynfam_inputs, override_settings


def extract_fam_tars(outputs_name, scratch=SCRATCH_DIR):
    """Extract every fam_meta/<OP-K>.tar produced for this nucleus in place,
    so the raw *_spin_*.dat files land where the analysis scripts expect them."""
    fam_meta = os.path.join(scratch, outputs_name, '000000', 'fam_soln', 'fam_meta')
    tars = sorted(glob.glob(os.path.join(fam_meta, '*.tar')))
    if not tars:
        print(f"  [warn] no tar files found in {fam_meta}")
        return
    for t in tars:
        with tarfile.open(t) as tf:
            tf.extractall(path=fam_meta)
        print(f"  extracted {os.path.relpath(t, scratch)}")


if __name__ == '__main__':
    # pynfamPaths.mkdirs() creates scratch/outputs with a single-level os.mkdir,
    # so SCRATCH_DIR itself must already exist (it's a new nested subfolder of
    # tests/, unlike the pre-existing top-level tests/ used by other scripts).
    if _RANK == 0:
        os.makedirs(SCRATCH_DIR, exist_ok=True)

    for Z, N in NUCLEI:
        A = Z + N
        label = f"{A}{ELEMENTS[Z]}"
        outputs_name = f"{label}_spindensity_test"

        if _RANK == 0:
            print(f"=== {label}  (Z={Z}, N={N})  ->  {outputs_name} ===")

        pynfam_inputs, override_settings = make_inputs(Z, N, outputs_name)
        pynfam_mpi_calc(pynfam_inputs, override_settings, check=False)

        if _RANK == 0:
            extract_fam_tars(outputs_name)
