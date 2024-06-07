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
    override_setts = {'vpair_t0_scaled':-0.7925, 'override_cs0': 129.297} #1.6, -1.7
    result = pynfam_fit_wrapper(input_data=filename, check=check, override_setts_fit=override_setts, scratch_dir='./fit_tests')
    print(result)
