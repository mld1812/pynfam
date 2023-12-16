# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
from future.builtins   import str
from future.builtins   import range
# -------------- Utilities -----------------
from copy import deepcopy
import os
import numpy as np
# ----------- Relative Imports -------------
from ..strength.phase_space import phaseSpace
from ..config import TMIN
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2021-06-10'


#===============================================================================#
#                             CLASS contourSetElem                              #
#===============================================================================#
class contourSetElem(object):
    """
    An instance of the class contourSetElem is an element of a set of contours.
    It contains a famContour object and addtional attributes and methods to
    identify it and adjust it (rotate, deform, shift bounds, etc.) as needed for
    calculations which use multiple contours per fam calculation (e.g. FT, ECFT, Q_eff).
    """
    id_tol = 1e-1 # Ignore differences less than 0.1
    FT_pole_label = u'FT_pole'

    def __init__(self, ctr, index, label1, label2=u""):
        self.ctr = ctr
        self.label1 = label1
        self.label2 = label2
        self.index = index
        self.all_strengths = None

        self.split_id = id(self)
        self.split_status = False

        self.is_ft_pole = False
        self.touching_ft_pole = False

    @property
    def subdir(self):
        return os.path.normpath(os.path.join(self.label1, self.label2))

    @property
    def ctr_id(self):
        """
        Round the bounds to a given precision
        We round the way for both bounds so the ids are continuous for splits
        We round up always to favor 1. E>0 and 2. Polynomial over Exponential PSECFT
        * NB: ctr_id only affects calculations with multiple different contours

        For invalid intervals, or intervals smaller than an order of magnitude
        less than the precision, return None to indicate this contour is invalid.
        """
        id_tol = contourSetElem.id_tol
        #bounds = (round(np.ceil(self.ctr.energy_min/id_tol)*id_tol, -int(np.log10(id_tol))),
        #          round(np.ceil(self.ctr.energy_max/id_tol)*id_tol, -int(np.log10(id_tol))))

        bounds = (round(self.ctr.energy_min, -int(np.log10(id_tol))),
                  round(self.ctr.energy_max, -int(np.log10(id_tol))))

        if self.is_ft_pole:
            return contourSetElem.FT_pole_label

        if bounds[0] >= bounds[1] or abs(bounds[0]-bounds[1]) < 1e-10 or any(np.isnan(bounds)):
            return None
        else:
            return bounds

    @property
    def ctr_valid(self):
        emin = self.ctr.energy_min
        emax = self.ctr.energy_max
        return not (emin >= emax or abs(emin-emax) < 1e-10 or any(np.isnan([emin, emax])))

    def split(self, label, val, lt=True):
        """
        label for what is being split on
        energy value at which to split the contour
        lt=True means the split is the one at lower energy, existing ctr is at higher energy
        """
        id_tol = contourSetElem.id_tol

        if self.ctr_id is None:
            return
        if not (self.ctr.energy_min <= val <= self.ctr.energy_max):
            return

        # If val is within tol of bounds, don't split but shrink
        #   1. This makes phase space fit work better...
        #      We'll get all of poly or exp, rather than one plus a small bit of the other
        #   2. This is essential for splitting FT contours at zero, otherwise the pole
        #      would be included inside the contour if Emax <= tol or Emin >= -tol
        if abs(val-self.ctr.energy_min) <= id_tol:
            self.ctr.updateSettings({u'energy_min': val})
            return
        elif abs(val-self.ctr.energy_max) <= id_tol:
            self.ctr.updateSettings({u'energy_max': val})
            return

        split = deepcopy(self)
        if lt:
            self.label2 = "spGT{:}".format(label)
            split.label2 = "spLT{:}".format(label)
            self.ctr.updateSettings({u'energy_min': val})
            split.ctr.updateSettings({u'energy_max': val})
        else:
            self.label2 = "spLT{:}".format(label)
            split.label2 = "spGT{:}".format(label)
            self.ctr.updateSettings({u'energy_max': val})
            split.ctr.updateSettings({u'energy_min': val})

        # Remember identity of the original object by overwriting id attribute
        split.split_id = self.split_id
        split.split_status = True
        self.split_status = True

        return split

    def rorate_ctr(self, rot):
        t0 = self.ctr._settings[u'theta_init']
        self.ctr.updateSettings({u'theta_init': t0-rot})

    def adjust_lb(self, val):
        """ If bounds are close to val, set them to exactly val """
        id_tol = contourSetElem.id_tol
        adjusted = False
        if abs(self.ctr.energy_min - val) < id_tol:
            self.ctr.updateSettings({u'energy_min': val})
            adjusted = True
        return adjusted

    def adjust_rb(self, val):
        id_tol = contourSetElem.id_tol
        adjusted = False
        if abs(self.ctr.energy_max - val) < id_tol:
            self.ctr.updateSettings({u'energy_max': val})
            adjusted = True
        return adjusted

    def set_bounds_to_id(self):
        if not isinstance(self.ctr_id, tuple):
            return
        self.ctr.updateSettings({u'energy_min': self.ctr_id[0]})
        self.ctr.updateSettings({u'energy_max': self.ctr_id[1]})

    def deform_ctr(self, temp):
        """
        The FT prefactor has poles all along imaginary axis, make sure edge of contour
        stays sufficiently far away from them by deforming to ellipse if necessary.
        hmax is max height such that contour is xtol away from 1st pole on imag axis.
        """
        if not isinstance(self.ctr_id, tuple):
            return

        # Use the rounded bounds to set H so uniqueness is still defined by ctr_id
        lb = self.ctr_id[0]
        rb = self.ctr_id[1]
        h  = self.ctr._settings[u'max_height']
        r0 = 0.5*(lb + rb)
        r = rb - r0

        # Location of 1st pole on imag axis used for distance away from imag axis
        hpole = 2*np.pi*temp
        xtol = np.sign(r0)*hpole

        # I'd rather get close to the pole than deform so much it becomes a strength fct
        # Don't deform ellipse any smaller than 1MeV
        hmin = 1

        # If contour doesn't cross x=xtol, theres no solution hmax to intersection w/ (xtol, hpole)
        if abs(r0 - xtol) >= r:
            # whole contour inside region, just make sure h < hpole
            if abs(r0) < abs(xtol):
                hmax = min(h, hpole)
            # whole contour outside region, no problems...
            else:
                hmax = h
        # Contour crosses x=xtol, make sure h at xtol < hpole
        else:
            hmax = min(h, hpole/np.sin(np.arccos((xtol - r0)/r)))

        hmax = max(hmax, hmin)
        self.ctr.updateSettings({u'max_height': hmax})

    def mkdirs(self, mgr, fam, label3=u""):
        if fam:
            p = mgr.paths.fam
            pm = mgr.paths.fam_m
        else:
            p = mgr.paths.beta
            pm = mgr.paths.beta_m
        if self.label1:
            mgr.paths.mkdirs([os.path.join(p, self.label1),
                              os.path.join(pm, self.label1)])
        if self.ctr_id is not None:
            if self.label2:
                mgr.paths.mkdirs([os.path.join(p, self.label1, self.label2),
                                  os.path.join(pm, self.label1, self.label2)])
            if label3:
                mgr.paths.mkdirs([os.path.join(p, self.label1, self.label2, label3),
                                  os.path.join(pm, self.label1, self.label2, label3)])

    def set_is_ft_pole(self):
        self.is_ft_pole = True
        self.label1 = contourSetElem.FT_pole_label

    def set_touching_ft_pole(self):
        self.touching_ft_pole = True

#-------------------------------------------------------------------------------
def get_contour_set(ctr_list):
    """ Get a list of unique contours"""
    ctr_set = list({c.ctr_id : deepcopy(c) for c in ctr_list if c.ctr_id is not None}.values())
    ctr_set.sort(key=lambda x: x.subdir)

    # If we are in a multi-contour mode with splitting, use the rounded bounds
    splits = [c for c in ctr_list if c.split_status]
    if len(ctr_set) > 1 and len(splits) > 0:
        for ctre in ctr_set:
            ctre.set_bounds_to_id()

    return ctr_set

#-------------------------------------------------------------------------------
def get_contour_split_set(ctr_list):
    """ Group contours based on split """
    id_set = list(set([c.split_id for c in ctr_list]))
    ctr_split_set = [[c for c in ctr_list if c.split_id==i] for i in id_set]
    return ctr_split_set

#-------------------------------------------------------------------------------
def multiple_psi_or_ctr(mgr, setts, hfb_gs, ctr_main):
    """
    PSI settings do not require re-running FAM or HFB. Allow any PSI override
    to be a list, in which case we re-run the beta calculation with every setting
    value in the list. psi_setting lists are treated the same as tuples are - every
    list provided must have the same length. The i'th element in the list outputs to
    a subdirectory in beta_soln named "psi_zfill(i,4)"

    log(pYe) and Q_eff are special cases where we might need to re-do the FAM calculation
    with a contour who's interval is adjusted for some or all of the elements in the
    psi input list. The FAM solutions are sent to subdirectories in fam_soln with the same
    name as those in beta_soln. This "internalized" workflow is messy, but I thinks it's the
    best option. I want to keep 1 master = 1 fam calculation as much as possible for memory
    reasons, only allowing multiples for closed contours in these special cases.
    """
    # Check for any lists as inputs. Doing it here instead of in workflow utils
    # lets us separate tuple functionality from list functionality
    psi_set = setts[u'psi']
    nr_psi = None
    for h in psi_set:
        if isinstance(psi_set[h], list):
            if len(psi_set[h]) == 0:
                raise RuntimeError(u"PSI input list with length zero detected")
            if nr_psi is None:
                nr_psi = len(psi_set[h])
            else:
                if len(psi_set[h]) != nr_psi:
                    raise RuntimeError(u"PSI input lists must have same length")
    if nr_psi is None: nr_psi = 1

    # Subdirectory labels
    labels = [u'']
    if nr_psi > 1:
        labels = [u'psi_{:}'.format('{:}'.format(i).zfill(3)) for i in range(nr_psi)]

    # Special cases where contours depend on psi_settings
    # So far: Q_eff and FTEC (FTEC can use Qeff for EQRPAmax...)
    # For now, Qeff is overriddent with log(pYe), rather than used in conjunction
    contours = None
    contours = multiple_qeff(labels, setts, hfb_gs, ctr_main)
    if hfb_gs.beta == u"c" and hfb_gs.ft_active:
        contours = ftec_multiple_densities(labels, setts, hfb_gs, ctr_main)

    # Handle all contours for this calc are the SAME
    # i.e. Same HFB, Same FAM, multiple BETA
    if contours is None:
        contours = [contourSetElem(deepcopy(ctr_main), i, lab) for i,lab in enumerate(labels)]

    # Handle finite temp calcs with closed contours
    if ctr_main.closed and hfb_gs.ft_active:
        contours = ft_contours(contours, setts, ctr_main, hfb_gs.temperature)

    return contours

#-------------------------------------------------------------------------------
def multiple_qeff(labels, setts, hfb_gs, ctr_main):
    """
    """

    # Initialize some variables
    pso = phaseSpace(hfb_gs.beta)
    pso.updateSettings(setts[u'psi'])
    eqrpamax = hfb_gs.soln_dict[u"EQRPA_max"]
    Q_hfb = hfb_gs.soln_dict[u"HFB_Qval"]
    Egs  = eqrpamax - Q_hfb

    # Force Q_eff to be iterable
    Q_effi = np.array(pso._settings[u'Q_eff'], copy=False, ndmin=1)
    Q_eff = pso._settings[u'Q_eff']
    Q_mode = int(pso._settings[u'Q_eff_mode'])

    # Every Q_eff uses same contour for Q_mode < 0
    if Q_eff is None or Q_mode <= 0:
        return

    # Adjust the max QRPA energy to use an effective Q value
    # NB: THIS SUPERCEDES ENERGY_MAX OVERRIDE...
    emin = ctr_main.energy_min
    emax_max = emin
    imax = 0
    contours = []
    for i, lab in enumerate(labels):
        # Mode 1: q_eff replaces q_hfb
        if Q_mode == 1:
            emax = (Q_effi[i]) + Egs
        # Mode 2: q_eff is a shift from q_hfb
        elif Q_mode == 2:
            emax = (Q_hfb + Q_effi[i]) + Egs

        if emax > emax_max:
            emax_max = emax
            imax = i

        ctr = deepcopy(ctr_main)
        ctr.updateSettings({u'energy_max': emax})
        contours.append(contourSetElem(ctr, i, lab))

    # Q_mode > 0, open contours use 1 ctr, but with largest emax
    if not ctr_main.closed:
        ctre = contours[imax]
        contours = [contourSetElem(deepcopy(ctre.ctr), i, lab) for i,lab in enumerate(labels)]

    return contours

#-------------------------------------------------------------------------------
def ftec_multiple_densities(labels, setts, hfb_gs, ctr_main):
    """ For finite temp EC calculations, we would like to compute the rates using different
    phase space corresponding to different stellar electron densities. The phase space becomes
    tricky for 2 main reasons:
    1. The phase space dies off rapidly to zero with QRPA energy depending on the density.
       A function with a large region of zero is poorly approximated by an analytic function.
    2. The character of the FT phase space function changes across the electron fermi energy,
       (which is a function of the density) from mostly polynomial to mostly exponential.
       It is difficult to reproduce this behavior with an analytic function.
    Both of these issues cause the quality of any analytic approximation to the phase space
    to be highly density-dependent. To combat this, we use multiple contours.

    For each density, we calculate on a contour on [0,Emax[rho]], where Emax[rho] is the maximum
    QRPA energy with non-zero phase space for a given density. If the fermi energy mu falls within
    this region, the contour is split into two contours [0,mu[rho]] and [mu[rho], Emax[rho]]. This
    procedure is carried out for each density, avoiding computing any repeat contours.
    """

    # Initialize some variables
    pso = phaseSpace(hfb_gs.beta)
    pso.updateSettings(setts[u'psi'])
    temper = hfb_gs.temperature
    (Z, N, A) = hfb_gs.nucleus
    eqrpamax = hfb_gs.soln_dict[u"EQRPA_max"]

    # Open contours - every density uses same contour
    if not ctr_main.closed:
        return

    emin = ctr_main.energy_min
    emax = ctr_main.energy_max
    emax_arr, mu_arr = pso.get_max_eqrpa_ecft(temper, Z, A, eqrpamax, emin, emax)
    if mu_arr.size != len(labels): raise RuntimeError

    # We related EQRPA to W0 in phase space. Wth= 1, or -W0 if W0<-1.
    # When W0<-1 = (mu/mec2+1) we have behavior change, or -W0 = (mu/mec2+1).
    # W0 = (Emax-EQRPA)/mec2 - 1, which gives EQRPAmu = Emax + mu
    emu_arr = eqrpamax + mu_arr

    contours = []
    for i, (lab, emu) in enumerate(zip(labels, emu_arr)):
        ctre1 = contourSetElem(deepcopy(ctr_main), i, lab)

        # ADJUST EMAX:
        #   Set emax to where phase space goes to zero.
        ctre1.ctr.updateSettings({u"energy_max": emax_arr[i]})

        ctre2 = ctre1.split(u"mu", emu)
        if ctre2 is not None:
            contours.append(ctre2)
        contours.append(ctre1)

    return contours

#-------------------------------------------------------------------------------
def ft_contours(contours, setts, ctr_main, temp):
    """
    Loop through all contours and split at zero if necessary. Also add the pole contour.

    For finite temp calcs, strength at positive and negative energies can contribute to the rate.
    However, the prefactor which makes this possible has a pole at zero.

    With closed contours, we need to avoid this pole at zero. The integration becomes numerically
    inaccurate if the contours are very close to but not intersecting the pole. If the contour is
    intersecting the pole, interestingly we pick up 1/2 the residue of the pole (calculus of partial residues)
    with no numerical issues, as long as the points are rotated such that none are exactly on the
    pole. My solution is to calculate 3 contours, (emin,0)+(0,emax)+(-tol,tol). The first two pick
    up the strength + 1/2 residue of the pole. The last one picks up only the residue of the pole,
    which we can then subtract out.

    ** The integration for a contour very close to but just outside the pole becomes better if we use trapezoidal
       rule with very specifically rotated points. The amount of rotation depends on the radius of the contour,
       the distance from the pole, and the number of points. But I havent been able to figure out what the magic
       number is analytically. It is NOT half the angle spacing as you might think.
    """

    pole_tol = 1e-6
    pref_tol = 1e-20

    # We override t0 here to ensure we get a sensible integration result with the FT pole
    # NB:
    # For npts= odd: point at t0 and t0+pi are on real axis - rotate by spacing/2 to avoid both
    # For npts=even: point at t0 is on real axis, point at t0+pi is maximally off
    #                rotating by spacing/2 puts t0 maximally off, but t0+pi on.
    #                However, any other rotation spoils the +/-i symmetry.
    t0shift = 0
    if not ctr_main.use_gauleg:
        npts = ctr_main.nr_points
        t0shift = np.pi/(npts-1) # 1/2*(2pi/nspaces)

    contours_ft = []
    for ctre1 in contours:
        emin = ctre1.ctr.energy_min
        emax = ctre1.ctr.energy_max
        lab2 = ctre1.label2
        lab1 = ctre1.label1
        pref = lab2

        if not ctre1.label2:
           label2 = u"0"
        else:
           label2 = u"0_{:}".format(ctre1.label2)

        # ADJUST EMIN AND THETA0:
        #   Strength at negative energies is important at FT, but the FT prefactor
        #   dies from -inf to 0 at more negative energies. We don't want to include portions
        #   where expontential in prefactor overflows and prefactor is exactly zero because
        #   the contour integration suffers, so use a cutoff. Some values for reference:
        #   FFN,tol=1e-20: [T=0.000861733, emin=-0.0397], [T=1, emin=-46.05], [T=8.61733, emin=-396.8]
        if temp <= TMIN:
            emin = max(emin, 0.0)
        else:
            emin = max(emin, temp*np.log(pref_tol/(1-pref_tol)))
        ctre1.ctr.updateSettings({u"energy_min": emin})

        # Don't use FT prefactor or E<0 below empirical cutoff temp where prefactor ~ unit step.
        if temp <= TMIN:
            contours_ft.append(ctre1)
            continue

        ctre2 = ctre1.split(label2, 0, lt=True)
        # For the split_lt0, rotate points by pi (so t0=pi becomes t0=0)
        if ctre2 is not None:
            ctre1.set_touching_ft_pole()
            ctre2.set_touching_ft_pole()
            ctre2.deform_ctr(temp)

            # Ensure there's never a point on (0,0) by rotating appropriately:
            # - Gauleg needs dense points near 0 (t0=pi for gt0, t0=0 for lt0)
            # - Equal spacing needs rotated by 1/2 spacing AND for even npts needs t0 same as Gauleg
            ctre1.ctr.updateSettings({u"theta_init": np.pi + t0shift})
            ctre2.ctr.updateSettings({u"theta_init": 0.0 + t0shift})

            contours_ft.append(ctre2)
        # If zero was not in the interval, identify if it is close to the bounds
        # This should be redundant with split behavior, but helps to flag bounds touching zero
        else:
            adj = ctre1.adjust_lb(0)
            if adj:
                ctre1.set_touching_ft_pole()
                ctre1.ctr.updateSettings({u"theta_init": np.pi + t0shift})
            else:
                adj = ctre1.adjust_rb(0)
                if adj:
                    ctre1.set_touching_ft_pole()
                    ctre1.ctr.updateSettings({u"theta_init": 0.0 + t0shift})
                # Contour not split at 0 or touching 0, but rotate same as split anyway
                # If it's close to 0 but not within tol, integration could still be poor
                # unless properly rotated (typically for small temps and densities)
                elif ctre1.ctr.energy_max < 0:
                    ctre1.ctr.updateSettings({u"theta_init": 0.0 + t0shift})
                else:
                    ctre1.ctr.updateSettings({u"theta_init": np.pi + t0shift})

        ctre1.deform_ctr(temp)
        contours_ft.append(ctre1)

    # Add the pole contour if we need to subtract it out
    do_ft = [c for c in contours_ft if c.touching_ft_pole]
    if do_ft:
        cpole = contourSetElem(deepcopy(ctr_main), None, u"FT_pole")
        cpole.set_is_ft_pole()
        cpole.ctr.updateSettings({'energy_min': -pole_tol, 'energy_max': pole_tol})
        contours_ft.append(cpole)

    return contours_ft
