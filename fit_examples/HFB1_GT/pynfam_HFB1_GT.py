#!/usr/bin/env python
import sys
sys.path.append("../..")
sys.path.append("..")
from pynfam import pynfam_fit_wrapper

# code for HFB1_GT run
if __name__ == '__main__':

    filename = 'pynfam_input_GT.csv'
    override_setts = {'g0p':1.6, 'g1p':0.0, 'h0p':0.0, 'vpair_t0_scaled':-1.7}
    pynfam_fit_wrapper(input_data=filename, override_setts_fit=override_setts, exes_dir="../../exes")
