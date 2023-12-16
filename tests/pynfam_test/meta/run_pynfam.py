#!/usr/bin/env python
from __future__ import unicode_literals
from pynfam import pynfam_mpi_calc

pynfam_inputs = {
 'directories': {
     'outputs' : 'pynfam_test',
     'exes'    : './exes',
     'scratch' : './tests'
     },

 'nr_parallel_calcs': None,

 'rerun_mode': 0,

 'hfb_mode': {
     'gs_def_scan'    : (0, ()),
     'dripline_mode'  : 0,
     'ignore_nonconv' : 0
     },

 'fam_mode': {
     'fam_contour': 'CONSTl',
     'beta_type'  : '-',
     'fam_ops'    : ('P', 0)

     }
}

override_settings = {
 'hfb' : {'proton_number'    : 36,
          'neutron_number'   : 80,
          'number_of_shells' : 6, 
          'number_gauss'     : 20, 
          'number_laguerre'  : 20, 
          'number_legendre'  : 40
          },

 'ctr' : {'nr_points' : 50},

 'fam' : {},
 #'two_body_current_mode' : 111
 'psi' : {}
}


if __name__ == '__main__':
    pynfam_mpi_calc(pynfam_inputs, override_settings, check=False)
