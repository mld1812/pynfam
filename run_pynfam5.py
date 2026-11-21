#!/usr/bin/env python
from __future__ import unicode_literals
from pynfam import pynfam_mpi_calc

pynfam_inputs = {
 'directories': {
     'outputs' : '230Hg_1-_test_new2bc', #  318Sg 106 212  
     'exes'    : './exes',
     'scratch' : './tests'
     },

 'nr_parallel_calcs': None,

 'rerun_mode': 0,

 'hfb_mode': {
     'gs_def_scan'    : (-2, (-0.4,-0.3,-0.2,-0.1, 0.0,0.1, 0.2, 0.3, 0.4)),
     'dripline_mode'  : 0,
     'ignore_nonconv' : 0
     },

 'fam_mode': {
     'fam_contour': 'CIRCLE',
     'beta_type'  : '-',
     'fam_ops'    : '1-'

     }
}


override_settings = {
 'hfb' : {'proton_number'    : 80,
          'neutron_number'   : 150,
          'number_of_shells' : 16, 
          #'number_gauss'     : 20,
          #'number_laguerre'  : 20,
          #'number_legendre'  : 40,
          'functional'       : 'SKOP',
          #'pairing_cutoff': 60.0,
          #'pairing_feature': 0.5,
          #"vpair_n": 0.0,
     	#"vpair_p": 0.0,
     	"number_iterations": 600,
     	"accuracy": 1e-08,
         #'proton_blocking': [-1,0,0,0,0],
        # 'neutron_blocking': [1,0,0,0,0]
          #'force_parity'     : False
          },

 'ctr' : {
        'nr_points': 60,        
        #'half_width': 0.01,
        #'de_hw_ratio': 0.25,
        #'hfb_emin_buff': 2.0
        #'energy_min': 1.75,
        #'energy_max': 1.8,
        },

 'fam' : {
    'use_fam_storage': 0,
        'interaction_name': 'SKOP',
        'force_j2_terms': False,
        'convergence_epsilon': 1e-07,
        'require_self_consistency': True,
        'require_gauge_invariance': True,
        'vpair_t1': None,
        #'override_cs0': 120.725, #Try this fit value, or 128.279 (1A)
        #'override_csr': 0.0,
        #'override_cds': 0.0,
        #'override_ct': None, # don't override tensor term
        #'override_cgs': None, # don't override tensor term
        #'override_cf': None, # don't override tensor term
        'max_iter': 300,
        'broyden_history_size': 200,

        #"vpair_t0": -361.74376445262425,
        #  "override_cs0": 118.46989547978822,
        'two_body_current_mode': 141221
 },

 'psi' : {
    'GA': -1.0,
    'UseNewShapeFactors': False
 }
}


if __name__ == '__main__':
    pynfam_mpi_calc(pynfam_inputs, override_settings, check=False)
