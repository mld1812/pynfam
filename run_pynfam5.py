#!/usr/bin/env python
from __future__ import unicode_literals
from pynfam import pynfam_mpi_calc
import numpy as np

#notes: openMP = 1. this is 2bc so we may want to compare with and without. 

pynfam_inputs = {
 'directories': {

     'outputs' : 'Ar48_rerun2Test',
     'exes'    : './exes',
     'scratch' : './tests'
     },
'nr_parallel_calcs': None,

 'rerun_mode': 2,

 'hfb_mode': {
     'gs_def_scan'    : (-2, (-0.2, 0.0, 0.2)),#(+1, ((0.24,)*400)),
     'dripline_mode'  : 0,
     'ignore_nonconv' : 0,
     },

 'fam_mode': {
     'fam_contour': 'CIRCLE',
     'beta_type'  : '-',
     'fam_ops'    : ('GT', 0), 

     }
}

override_settings = {
 'hfb' : {'proton_number'    : 18, 
          'neutron_number'   : 30,
          'number_of_shells' : 8,
          'number_gauss'     : 40,
          'number_laguerre'  : 40,
          'number_legendre'  : 80,
	        #'functional'       : 'HFB1',
          #'force_parity'     : False,
          },


 'ctr' : {#'energy_min' : 2.4,
          #'energy_max' : 22.4,
          #'half_width' : 0.2,
          #'nr_points'  : 201,

     },

 'fam' : {
     'two_body_current_mode' : 111100,
     #'interaction_name' : 'SKM*',
          #'vpair_t0'         : 0.0,
          #'convergence_epsilon' : 1e-2,
          #'max_iter'         : 100,
          #'require_gauge_invariance': False,
          #'override_cs0'        : '',
          #'override_csr'        : '',
          #'override_cds'        : '',
          #'override_ct'         : '',
          #'override_cgs'        : '',
          #'override_cf'         : '',
 
     },

 'psi' : {
	#'GA': -1.0
}

}


if __name__ == '__main__':
    pynfam_mpi_calc(pynfam_inputs, override_settings, check=False)
