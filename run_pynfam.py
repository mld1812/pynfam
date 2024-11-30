#!/usr/bin/env python
from __future__ import unicode_literals
from pynfam import pynfam_mpi_calc

pynfam_inputs = {
 'directories': {
     'outputs' : '0-_operator_16shell_1111_test4',
     'exes'    : './exes',
     'scratch' : './tests'
     },

 'nr_parallel_calcs': None,

 'rerun_mode': 0,

 'hfb_mode': {
     'gs_def_scan'    :  (0, ()), #(-2, (-0.2, 0, 0.2)), #(0, ()), #+1, (0.14,)
     'dripline_mode'  : 0,
     'ignore_nonconv' : 0,
     #'retry_nonconv': True,
     },

 'fam_mode': {
     'fam_contour': 'CIRCLE',
     'beta_type'  : '-',
     'fam_ops'    : ('RS0i', 0) #('P', 0)
     }
}

override_settings = {
 'hfb' : {'proton_number'    : 16, #64, 98 - Gd162
          'neutron_number'   : 24,
          'number_of_shells' : 16, 
          #'number_gauss'     : 20,
          #'number_laguerre'  : 20,
          #'number_legendre'  : 40,
          'functional'       : 'SKOP',
          'pairing_cutoff': 60.0,
          'pairing_feature': 0.5,
          #"vpair_n": 0.0,
     	#"vpair_p": 0.0,
     	"number_iterations": 600,
     	"accuracy": 1e-08,
          #'force_parity'     : False
          },

 'ctr' : {
	    #'energy_min': 5.287984000000005,
	    #'energy_max': 39.287984,
	    #'half_width': 0.5, #0.5
	    #'de_hw_ratio': 0.8
	    'nr_points' : 6
	   },

 'fam' : {
 	 'use_fam_storage': 0,
 	 'interaction_name' : 'SKOP',
 	 'force_j2_terms': False,
 	 'convergence_epsilon': 1e-07,
 	 'require_self_consistency': True,
 	 'require_gauge_invariance': False,
 	 #'vpair_t0': -264.216261162,
 	 'vpair_t1': None,
 	 #'override_cs0': 137.00212941970156,
 	 'override_csr': 0.0,
 	 'override_cds': 0.0,
 	 'override_ct' : 0.0,
 	 'override_cgs' : 0.0,
 	 'override_cf': 0.0,
 	 'max_iter': 300,
 	 'broyden_history_size': 200,
 
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
 
 	 #'two_body_current_mode': 111022, #digit 4: GT, digit 5: P, digit 6: PS0; 1 = 2bc, 2 = DME for P/PS0. 2nd digit controls GT: use 1 for full pnfam, 4 for DME
		#'two_body_current_lecs': [-9.5886, 0, 0],
 	  #'max_iter' : 2
          #'two_body_current_usep': False,
          #'two_body_current_lecs': [-3.1962, 3.1962, 0]
          },

 'psi' : {
 	#"GA": -1.0
 	#            (u'Q_eff'      , None),
        #    (u'Q_eff_mode' , 0),
        #'Q_eff': 7.0,
        #'Q_eff_mode': 1
 }

}


if __name__ == '__main__':
    pynfam_mpi_calc(pynfam_inputs, override_settings, check=False)
