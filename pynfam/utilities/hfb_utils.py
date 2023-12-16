# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
# ---------------- Utilities ---------------
from shutil import copy2
from shutil import rmtree
from copy import deepcopy
import os
import tarfile
import numpy as np
# ----------- Relative Imports -------------
from ..outputs.ls_logger import hfbLogger
from ..outputs.hfbtho_parser import hfbthoParser
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-08-28'


# ------------------------------------------------------------------------------
def hfbEvenDefList(hfb_main, kickoff, scan_list):
    """
    Generate a list of hfbthoRun objects for even-even nuclei with
    different deformations from a template hfbthoRun object.
    (Blocking is turned off in the resulting objects).

    Args:
        hfb_main (hfbthoRun): The hfbthoRun object template.
        kickoff (int): Constraint mode for HFBTHO (+/-2 is allowed as well)
        scan_list (tuple): List of beta2 deformations.

    Returns:
        list of hfbthoRun
    """

    # Adjust for kickoff=+/-2 supplied directly to hfb override settings
    t = deepcopy(hfb_main)
    t.setNmlParam({u'lambda_active':list(np.sign(t.activeConstraints))})

    # Don't do def scan if kickoff=0 or no defs in scan_list
    if not scan_list or not kickoff:
        tasks = [t]
    else:
        tasks = []
        # If we entered only 1 def and forgot the comma to make it a tuple...
        try:
            iter(scan_list)
        except TypeError:
            scan_list = (scan_list,)
        for basis_def in scan_list:
            hfb_obj = deepcopy(t)
            # Set basis def and initial woods saxon def
            hfb_obj.setBeta2s(basis_def)

            # Override any Q2 constraint for scan (and adjust for scan kickoff=+/-2)
            Q2 = hfb_obj.calculatedQ2() # NB: This uses odd mass for even core
            hfb_obj.setConstraint(1, Q2, np.sign(kickoff))
            # Add the kickoff to the soln dict so it appears in logfile
            tasks.append(hfb_obj)

    # Turn off blocking and temperature and set labels
    for i, nuc in enumerate(tasks):
        nuc.setBlocking(0, 0, 0)
        nuc.setTemperature(0.0)
        label_prefix = hfb_main.paths.rp(hfb_main.paths.hfb_m)
        if hfb_main.blocking[0][0] or hfb_main.blocking[1][0]:
            nuc.label = os.path.join(label_prefix, str(i).zfill(4), u'core')
        elif hfb_main.ft_active:
            nuc.label = os.path.join(label_prefix, str(i).zfill(4), u'zeroT')
        else:
            nuc.label = os.path.join(label_prefix, str(i).zfill(4))

    return tasks

#--------------------------------------------------------------
def hfbEvenRerunList(finished_nucs):
    """
    Generate a list of non-converged and converged hfbthoRun objects from a
    list of completed hfbthoRun objects. The basis deformation of the
    non-converged solutions is adjusted to the most recent value in the
    calculation, which can help reach convergence with a 2nd run.

    Args:
        finished_nucs (list): List of completed hfbthoRuns

    Returns:
        list of hfbthoRun: adjusted non-converged hfbthoRun objects
        list of hfbthoRun: converged hfbthoRun objects
    """

    # Check against 'No' so 'Error' does not get rerun
    nc  = [n for n in finished_nucs if n.soln_dict[u'Conv'] == u'No']
    fin = [n for n in finished_nucs if n.soln_dict[u'Conv'] != u'No']
    for i, nuc in enumerate(nc):
        new_basis = round(nuc.soln_dict[u'Def'],3)
        cur_basis = nuc.nml[u'HFBTHO_GENERAL'][u'basis_deformation']
        # If the basis is too similar, don't bother re-running...
        if abs(new_basis - cur_basis) < 0.001:
            fin.append(nc.pop(i)) # NB: remove nuc from nc!
        else:
            nuc.setRestartFile(False)
            nuc.setBeta2s(new_basis)

    return nc, fin

#--------------------------------------------------------------
def hfbOddList(even_objs, hfb_main, kickoff):
    """
    Generate a list of hfbthoRun objects for odd nuclei from a list of
    completed even-even solutions. All blocking candidates are considered,
    odd-odd nuclei consider every combination of proton and neutron candidates.

    Args:
        even_objs (list): List of completed (even) hfbthoRuns.
        hfb_main (hfbthoRun) : The original template object.
        kickoff (int): The constraint type for HFBTHO (+/- 2 is allowed).

    Returns:
        list of hfbthoRun
    """

    # Adjust for kickoff=+/-2 supplied directly to hfb override settings
    # NB: kickoff +/-1 only applies to even, +/-2 applies to both even and odd
    ac = list(np.sign([0 if abs(c)==1 else c for c in hfb_main.activeConstraints]))

    og_blocking = hfb_main.blocking
    tasks = []
    for core in even_objs:
        blocking_list = []
        # Get candidates and adjust omega for correct sign
        soln = hfbthoParser(core.output)
        candidates = soln.getBlockingCands()
        for i, cand_list in enumerate(candidates):
            for cand in cand_list:
                # 2*omega same sign as original nml for correct blocking
                cand[0] = cand[0]*int(np.sign(og_blocking[i][0]))

        # So double for-loop works for all cases, always give a list
        #  -- Single scan blocking: nucleon that is off is a 1 element list of [None]*5
        #  -- No scan blocking: 1 element list of og_blocking
        n_candidates, p_candidates = [5*[None]], [5*[None]]
        if og_blocking[0][0]:
            if not og_blocking[0][1]:
                n_candidates = candidates[0]
            else:
                n_candidates = [og_blocking[0]]
        if og_blocking[1][0]:
            if not og_blocking[1][1]:
                p_candidates = candidates[1]
            else:
                p_candidates = [og_blocking[1]]

        # Generate all blocking run objects
        count = 0
        for n_cand in n_candidates:
            for p_cand in p_candidates:
                blocked_obj = deepcopy(core)
                blocked_obj.pmt_inputs.append(os.path.join(core.rundir, core.file_bin))
                blocked_obj.setBlocking(n_cand, p_cand)
                blocked_obj.setRestartFile(True)
                # If no scan, set constraint from hfb_override setting
                # If scan kickoff==+/-1, remove constraint
                # If scan kickoff==+/-2, use same constraint as for core
                if not kickoff:
                    blocked_obj.setNmlParam({u'lambda_active':ac})
                elif abs(kickoff)==1:
                    blocked_obj.setConstraint(1, 0.0, 0)
                # Remove the "core" label replace it with "count"
                blocked_obj.core_rundir = core.rundir
                blocked_obj.label = os.path.join(os.path.dirname(core.label),
                        str(count).zfill(4))
                blocking_list.append(blocked_obj)
                count += 1

        tasks += blocking_list

    return tasks

#--------------------------------------------------------------
def hfbFTList(even_objs, hfb_main):
    og_temp = hfb_main.temperature
    tasks = []
    for count, core in enumerate(even_objs):
        ft_obj = deepcopy(core)
        ft_obj.setTemperature(og_temp)
        ft_obj.core_rundir = core.rundir
        ft_obj.label = os.path.join(os.path.dirname(core.label),
                        str(0).zfill(4))
        # Start up from 0T soln for efficiency
        ft_obj.pmt_inputs.append(os.path.join(core.rundir, core.file_bin))
        ft_obj.setRestartFile(True)
        tasks.append(ft_obj)
    return tasks

#--------------------------------------------------------------
def convString(conv_list):
    """
    Determine convergence string given a list of individual strings.

    String options are:

          * 'Yes'  : all convs = Yes
          * 'No'   : all convs = No
          * 'ATTN!': some convs=Yes, some convs=No
          * 'Error': There was a problem, or no convs=Yes or No

    Args:
        conv_list (list of str): The individual convergence strings.

    Returns:
        str
    """

    try:
        conv_yes = [conv_str for conv_str in conv_list if conv_str == u'Yes']
        conv_no  = [conv_str for conv_str in conv_list if conv_str != u'Yes']
        if conv_yes:
            conv = u'Yes'
            if conv_no:
                conv = u'ATTN!'
        elif conv_no:
            conv = u'No'
        else:
            raise RuntimeError
    except Exception:
        conv = u'Error'
    return conv

#--------------------------------------------------------------
def hfbConvCheck(hfb_list, ignore_nc):
    """
    Wrapper for hfbConvCheck to generate list of converged hfbthoRun
    objects and determine converged error status based on the pynfam
    input 'ignore_nc'.

    Args:
        hfb_list (list): List of completed hfbthoRuns.
        ignore_nc (int): The input file parameter.

    Returns:
        list: Properly converged solutions according to ignore_nc.
        bool: Any non converged solutions according to ignore_nc.

    Raises:
        ValueError
    """

    # Checking evens or odds?
    evens = [n for n in hfb_list if not (n.blocking[0][0] or n.blocking[1][0])]
    odds  = [n for n in hfb_list if (n.blocking[0][0] or n.blocking[1][0])]
    if evens and odds:
        msg = u"Must provide only even or only odd nuclei to pynfamConvCheck."
        raise ValueError(msg)

    # Use NC solutions
    if ignore_nc == 3:
        return hfb_list, False

    # Check entire list or only converged, depending on mode
    check_all = (evens and ignore_nc<2) or (odds and ignore_nc<1)
    if check_all:
        obj = hfb_list
    else:
        obj = [n for n in hfb_list if n.soln_dict[u'Conv'] == u'Yes']

    # Check convergence of the paired down list
    if not obj:
        conv_err = True
    else:
        conv_list= [hfb.soln_dict[u'Conv'] for hfb in obj]
        conv = convString(conv_list)
        conv_err = conv != u'Yes'

    return obj, conv_err

#-------------------------------------------------------------------------------
def hfbGroundState(hfb_list):
    """
    Determine the hfb ground state from a list hfbthoRun objects or
    hfbthoRun.soln_dict's.

    The ground state is determined by the lowest binding energy solution,
    but the deformation is also considered. We can get super-deformed
    meta-stable states with extremely low bindings energies, which are not
    true ground state solutions. We therefore require deformation beta2 be
    less than 0.6. Additionally, for odd nuclei if the lowest energy state
    is a super-deformed meta-stable state, we neglect the entire group of
    solutions for that parent nucleus. This ensures the chosen ground state
    is a true local minimum in the deformation+blocking_candidate space.

    Args:
        hfb_list (list): List of completed hfbthoRuns.

    Returns:
        hbfthoRun, dict, None: Return is same type as elements of input list
    """

    # Default return for empty input
    if not hfb_list:
        return None

    # Group by parent nucleus
    obj_labels = []
    for o in hfb_list:
        up1  = os.path.basename(os.path.dirname(o.label))
        base = os.path.basename(o.label)
        lab  = base
        if up1 != o.paths.subdir_hfb_m:
            lab = up1
        obj_labels.append({u'obj':o,u'lab':lab})
    labels = list(set([d[u'lab'] for d in obj_labels]))
    groups = [[d[u'obj'] for d in obj_labels if d[u'lab']  == l] for l in labels]

    # Deformation check
    max_def = 0.6
    good = []
    for g in groups:
        o_defs = [o for o in g if o.soln_dict[u'Def'] <= max_def]
        if not o_defs:
            min_def = None
        else:
            min_def = min(o_defs, key=lambda x:x.soln_dict[u'Energy'])
        min_all = min(g, key=lambda x:x.soln_dict[u'Energy'])
        if min_def == min_all:
            good.append(min_def)

    if not good:
        hfb_gs = None
    else:
        hfb_gs_og = min(good, key=lambda x:x.soln_dict[u'Energy'])
        hfb_gs = deepcopy(hfb_gs_og)

    return hfb_gs
