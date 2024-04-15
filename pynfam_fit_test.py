#!/usr/bin/env python
from pynfam import pynfam_fit_wrapper
import sys
import pickle

# code for test
if __name__ == '__main__':

    # input_data and check can be give by command line arguments
    # input_data should have ".csv" in it
    # check should be 0 or 1 or other ints
    # return is pickled for further check
    filename = 'pynfam_fit_input.csv'
    check = False
    if len(sys.argv)>1:
        if sys.argv[1][-4:] == '.csv':
            filename = sys.argv[1]
            if len(sys.argv)>2:
                check = bool(int(sys.argv[2]))
        else:
            check = bool(int(sys.argv[1]))
            
    # nr_points is set to a small number below for a quick test
    # g0p, g1p, h0p, vpair_t0 are fitting variables
    override_setts = {'g0p':1.6, 'g1p':0.0, 'h0p':0.0, 'vpair_t0_scaled':-1.7, 'nr_points':60}
    result = pynfam_fit_wrapper(input_data=filename, check=check, override_setts_fit=override_setts, scratch_dir='./fit_tests2')
    print(result)
