#!/usr/bin/env python
from __future__ import unicode_literals
from pynfam import pynfam_mpi_calc

pynfam_inputs = {
 'directories': {
     'outputs' : 'S40_RS0_2bc_onepoint',
     'exes'    : './exes',
     'scratch' : './tests'
     },

 'nr_parallel_calcs': None,

 'rerun_mode': 0,

 'hfb_mode': {
     'gs_def_scan'    : (0, ()), #+1, (0.14,)
     'dripline_mode'  : 0,
     'ignore_nonconv' : 0
     },

 'fam_mode': {
     'fam_contour': 'CONSTL',
     'beta_type'  : '-',
     'fam_ops'    : ('RS0', 0) #('P', 0)
     }
}

override_settings = {
 'hfb' : {'proton_number'    : 16, #64, 98 - Gd162
          'neutron_number'   : 24,
          'number_of_shells' : 6, 
          'number_gauss'     : 20,
          'number_laguerre'  : 20,
          'number_legendre'  : 40,
          #'functional'       : 'SKM*',
          #'force_parity'     : False
          },

 'ctr' : {
	    'energy_min': 2.4,
	    'energy_max': 22.4,
	    'half_width': 0.2, #0.5
	    'nr_points' : 2
	   },

 'fam' : {# 'interaction_name' : 'SKM*',
          #'vpair_t0'         : 0.0,
          #'convergence_epsilon' : 1e-2,
          'max_iter'         : 2,
          #'require_gauge_invariance': False,
          #'override_cs0'        : '',
          #'override_csr'        : '',
          #'override_cds'        : '',
          #'override_ct'         : '',
          #'override_cgs'        : '',
          #'override_cf'         : '',
 
 	 'two_body_current_mode': 111100, #digit 4: GT, digit 5: P, digit 6: PS0; 1 = 2bc, 2 = DME for P/PS0. 2nd digit controls GT: use 1 for full pnfam, 4 for DME
		
 	  #'max_iter' : 2
          #'two_body_current_usep': False,
          #'two_body_current_lecs': [-3.1962, 3.1962, 0]
          },

 'psi' : {}
}


if __name__ == '__main__':
    pynfam_mpi_calc(pynfam_inputs, override_settings, check=False)
