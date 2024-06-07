#!/usr/bin/env python
#Test spin-dipole resonance. 
from __future__ import unicode_literals
from pynfam import pynfam_mpi_calc

pynfam_inputs = {
 'directories': {
     'outputs' : 'pynfam_fit_wrapper_test',
     'exes'    : './exes',
     'scratch' : './tests'
     },

 'nr_parallel_calcs': None,

 'rerun_mode': 0,

 'hfb_mode': {
     'gs_def_scan'    :  (-2, (0.0,)), #(-2, (-0.2, 0, 0.2)), #(0, ()), #+1, (0.14,)
     'dripline_mode'  : 0,
     'ignore_nonconv' : 2,
     'retry_nonconv': True,
     },

 'fam_mode': {
     'fam_contour': 'CONSTR',
     'beta_type'  : '-',
     'fam_ops'    : 'SPINDIPOLE'#('P', 0)
     }
}

override_settings = {
 'hfb' : {'proton_number'    : 40, #64, 98 - Gd162
          'neutron_number'   : 50,
          'number_of_shells' : 16, 
          #'number_gauss'     : 20,
          #'number_laguerre'  : 20,
          #'number_legendre'  : 40,
          'functional'       : 'SKOP',
          'pairing_cutoff': 60.0,
          'pairing_feature': 0.5,
          "vpair_n": 0.0,
     	"vpair_p": 0.0,
     	"number_iterations": 600,
     	"accuracy": 1e-08,
          #'force_parity'     : False
          },

 'ctr' : {
	    'energy_min': 5.287984000000005,
	    'energy_max': 39.287984,
	    'half_width': 0.5, #0.5
	    'de_hw_ratio': 0.8
	    #'nr_points' : 2
	   },

 'fam' : {
 	 'use_fam_storage': 0,
 	 'interaction_name' : 'SKOP',
 	 'force_j2_terms': False,
 	 'convergence_epsilon': 1e-07,
 	 'require_self_consistency': True,
 	 'require_gauge_invariance': True,
 	 'vpair_t0': -264.216261162,
 	 'vpair_t1': None,
 	 'override_cs0': 153.31150462285333,
 	 'override_csr': 0.0,
 	 'override_cds': 0.0,
 	 'override_ct' : None,
 	 'override_cgs' : None,
 	 'override_cf': None,
     #'override_cj': None,
     #'override_csdj': None,
 	 'max_iter': 300,
 	 'broyden_history_size': 200
 
          #'vpair_t0'         : 0.0,
          #'convergence_epsilon' : 1e-2,
         # 'max_iter'         : 2,
          #'require_gauge_invariance': False,
          #'override_cs0'        : '',
          #'override_csr'        : '',
          #'override_cds'        : '',
          #'override_ct'         : '',
          #'override_cgs'        : '',
          #'override_cf'         : '',
 
 	 #'two_body_current_mode': 111020, #digit 4: GT, digit 5: P, digit 6: PS0; 1 = 2bc, 2 = DME for P/PS0. 2nd digit controls GT: use 1 for full pnfam, 4 for DME
		
 	  #'max_iter' : 2
          #'two_body_current_usep': False,
          #'two_body_current_lecs': [-3.1962, 3.1962, 0]
          },

 'psi' : {
 	"GA": -1.0
 }

}


if __name__ == '__main__':
    pynfam_mpi_calc(pynfam_inputs, override_settings, check=False)