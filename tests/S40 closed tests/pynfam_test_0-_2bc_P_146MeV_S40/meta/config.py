"""
This module contains the default values for hfbtho, pnfam,
contour, and phase space objects, as well as some commonly
used parameters.
"""
# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
# -------------- Utilities -----------------
from collections import OrderedDict as OD
import numpy as np
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2020-05-13'

DEFAULTS = {

    u'hfb': OD([
        (u'HFBTHO_GENERAL', OD([
            (u'number_of_shells'    , 16),
            (u'oscillator_length'   , -1.0),
            (u'basis_deformation'   , 0.0),
            (u'proton_number'       , 24),
            (u'neutron_number'      , 26),
            (u'type_of_calculation' , 1)
            ])),
        (u'HFBTHO_INITIAL', OD([
            (u'beta2_deformation'   , 0.0),
            (u'beta3_deformation'   , 0.0),
            (u'beta4_deformation'   , 0.0)
            ])),
        (u'HFBTHO_ITERATIONS', OD([
            (u'number_iterations'   , 300),
            (u'accuracy'            , 1e-06),
            (u'restart_file'        , 1)
            ])),
        (u'HFBTHO_FUNCTIONAL', OD([
            (u'functional'          , u'SKOP'),
            (u'add_initial_pairing' , False),
            (u'type_of_coulomb'     , 2)
            ])),
        (u'HFBTHO_PAIRING', OD([
            (u'user_pairing'        , False),
            (u'vpair_n'             , -250.0),
            (u'vpair_p'             , -250.0),
            (u'pairing_cutoff'      , 60.0),
            (u'pairing_feature'     , 0.5)
            ])),
        (u'HFBTHO_CONSTRAINTS', OD([
            (u'lambda_values'       , [1, 2, 3, 4, 5, 6, 7, 8]),
            (u'lambda_active'       , [0, 0, 0, 0, 0, 0, 0, 0]),
            (u'expectation_values'  , [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
            ])),
        (u'HFBTHO_BLOCKING', OD([
            (u'proton_blocking'     , [0, 0, 0, 0, 0]),
            (u'neutron_blocking'    , [0, 0, 0, 0, 0])
            ])),
        (u'HFBTHO_PROJECTION', OD([
            (u'switch_to_tho'       , 0),
            (u'projection_is_on'    , 0),
            (u'gauge_points'        , 1),
            (u'delta_z'             , 0),
            (u'delta_n'             , 0)
            ])),
        (u'HFBTHO_TEMPERATURE', OD([
            (u'set_temperature'     , False),
            (u'temperature'         , 0.0)
            ])),
        (u'HFBTHO_FEATURES', OD([
            (u'collective_inertia'  , False),
            (u'fission_fragments'   , False),
            (u'pairing_regularization', False),
            (u'automatic_basis'     , False),
            (u'localization_functions', False)
            ])),
        (u'HFBTHO_TDDFT', OD([
            (u'filter'              , False),
            (u'fragment_properties' , False),
            (u'real_z'              , 49.7),
            (u'real_n'              , 71.2)
            ])),
        (u'HFBTHO_NECK', OD([
            (u'set_neck_constrain'  , False),
            (u'neck_value'          , 13.0)
            ])),
        (u'HFBTHO_DEBUG', OD([
            (u'number_gauss'        , 40),
            (u'number_laguerre'     , 40),
            (u'number_legendre'     , 80),
            (u'compatibility_hfodd' , False),
            (u'number_states'       , 500),
            (u'force_parity'        , True),
            (u'print_time'          , 0)
            ])),
        (u'HFBTHO_RESTORATION', OD([
            (u'PNP_is_on'                  , 0),
            (u'number_of_gauge_points'     , 9),
            (u'delta_neutrons'             , 6),
            (u'delta_protons'              , 6),
            (u'AMP_is_on'                  , 0),
            (u'number_of_rotational_angles', 27),
            (u'maximal_angular_momentum'   , 10)
            ]))
        ]),

    u'fam': OD([
        (u'GENERAL', OD([
            (u'use_fam_storage'     , 0)
            ])),
        (u'EXT_FIELD', OD([
            (u'compute_crossterms'  , True),
            (u'two_body_current_mode', 0),
            (u'two_body_current_usep', False),
            (u'two_body_current_lecs', [-3.1962, 3.1962, 0])
            ])),
        (u'INTERACTION', OD([
            (u'interaction_name'    , u'SKOP'),
            (u'require_self_consistency', True),
            (u'require_gauge_invariance', True),
            (u'force_j2_terms'      , False),
            (u'vpair_t0'            , -346.352),
            (u'vpair_t1'            , None),
            (u'override_cs0'        , 128.279),
            (u'override_csr'        , 0.0),
            (u'override_cds'        , 0.0),
            (u'override_ct'         , None),
            (u'override_cgs'        , None),
            (u'override_cf'         , None)
            ])),
        (u'SOLVER', OD([
            (u'max_iter'            , 300),
            (u'convergence_epsilon' , 1e-07),
            (u'broyden_history_size', 50),
            (u'energy_shift_prot'   , 0.0),
            (u'energy_shift_neut'   , 0.0),
            (u'quench_residual_int' , 1.0)
            ]))
        ]),

    u'ctr': OD([
        (u'INTERVAL', OD([
            (u'energy_min'   , 0.0),
            (u'energy_max'   , 6.0),
            (u'hfb_emin_buff', 0.0)
            ])),
        (u'CIRCLE', OD([
            (u'nr_points'      , 60),
            (u'use_gauleg_ctr' , True),
            (u'beta_quadrature', u'GAUSS'),
            (u'shift_imag'     , 0.0),
            (u'theta_init'     , np.pi),
            (u'max_height'     , 30)
            ])),
        (u'CONSTL', OD([
            (u'nr_points'       , 60),
            (u'half_width'      , 0.1),
            (u'beta_quadrature' , u'TRAP')
            ])),
        (u'CONSTR', OD([
            (u'half_width'      , 0.1),
            (u'de_hw_ratio'     , 1.0),
            (u'beta_quadrature' , u'TRAP')
            ])),
        (u'EXP', OD([
            (u'hw_min'             , 0.01),
            (u'hw_max'             , 0.4),
            (u'p_percent_interval' , -0.25),
            (u'de_hw_ratio'        , 1.0),
            (u'beta_quadrature'    , u'TRAP')
            ])),
        (u'MONOMIAL', OD([
            (u'hw_min'             , 0.01),
            (u'hw_max'             , 0.4),
            (u'power'              , 1.0),
            (u'de_hw_ratio'        , 1.0),
            (u'beta_quadrature'    , u'TRAP')
            ])),
        (u'FERMIS', OD([
            (u'hw_min'             , 0.01),
            (u'hw_max'             , 0.40),
            (u'u_percent_interval' , 0.50),
            (u't_percent_interval' , 0.10),
            (u'de_hw_ratio'        , 1.0),
            (u'beta_quadrature'    , u'TRAP')
            ])),
        (u'FERMIA', OD([
            (u'nr_points_min'      , 10),
            (u'nr_points_max'      , 70),
            (u'hw_min'             , 0.01),
            (u'hw_max'             , 0.40),
            (u'u_percent_interval' , 0.50),
            (u't_percent_interval' , 0.10),
            (u'de_hw_ratio'        , 1.0),
            (u'beta_quadrature'    , u'TRAP')
            ]))
        ]),

    u'psi': OD([
        (u'PSI', OD([
            (u'psi_approx' , 'RATINT'),
            (u'psi_glpts'  , 15),
            (u'screening'  , None),
            (u'Q_eff'      , None),
            (u'Q_eff_mode' , 0),
            (u'GA'         , -1.0),
            (u'GV'         , 1.0),
            (u'log(pYe)'   , 9)
            ])),
        (u'RATINT', OD([
            (u'ratint_pts', 20)
            ])),
        (u'POLYFIT', OD([
            (u'polynomial_order'   , 6),
            (u'polynomial_fit_pts' , 200)
            ]))
        ])
}

# Constants
ALPHA    = 7.2973525698e-03 # fine-structure const.
HBARC    = 197.3269718      # hbar c [MeV fm]
HBAR_MEC = 386.15926800     # (hbar c)/(me c^2) [fm]
MEC2     = 0.510998928      # me c^2 [MeV]
KAPPA    = 6147.0           # decay const.
DMASS_NH = 0.78227          # n-H mass difference [MeV]
R0       = 1.2              # nuclear radius [fm]
MN       = 939.0            # nucleon mass [MeV]
NA       = 6.022140858e23   # Avogadros number [1/mol]
KB       = 8.6173304e-11    # Boltzmann constant [MeV/K]
TMIN     = 1.0*1e9*KB       # Temperature cutoff [MeV] below which FT prefactor is assumed a unit step

EPSILON  = np.finfo(np.float_).eps # machine precision for double
