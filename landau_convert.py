# convert between the experimental Landau parameter (pi+rho+g' model) and the one used in HFODD / PNFAM
# Author: Tong Li, Mich State Univ, 2020-
import numpy as np
from pathlib import Path

inverse_n0_exp = 392 # 1/N'_0, used in exp papers

def landau_from_exp(g_exp, rho_nm, smass_nm, hbzero):
    # g_exp is the experimental Landau parameter, 
    # see Sec. 3.4 of Progress in Particle and Nuclear Physics 56 (2006) 446
    # or Physics Reports 155, No. 5 (1987) 263
    # rho_nm is the saturation density. 
    # smass_nm is m / m*, where m* is the effective mass. 
    # rho_nm, smass_nm and hbzero are given in HFBTHO output. 
    # return is the Landau parameter following the convention of HFODD / PNFAM
 
    n0 = (3/2 * rho_nm)**(1/3) / (np.pi**(4/3) * hbzero * smass_nm)  # N_0 used in HFODD / PNFAM
    return g_exp * inverse_n0_exp * n0

def landau_to_exp(g, rho_nm, smass_nm, hbzero):
    # g is the Landau parameter used in HFODD / PNFAM, 
    # rho_nm is the saturation density. 
    # smass_nm is m / m*, where m* is the effective mass. 
    # rho_nm, smass_nm and hbzero are given in HFBTHO output. 
    # return is the Landau parameter following the convention of experimental papers. 
    
    inverse_n0 = np.pi**(4/3) * hbzero * smass_nm / ((3/2 * rho_nm)**(1/3))  # 1/N_0 used in HFODD / PNFAM
    return g * inverse_n0 / inverse_n0_exp

def read_hfb(output_file):
    # read functional info from HFB calcs
    functional_str = {'RHO_NM=', 'SMASS_NM=', 'hbzero='}
    functional_info = {}
    for line in output_file:
        found = set()
        for try_str in functional_str:
            pos = line.find(try_str)
            if pos >= 0:
                found.add(try_str)
                list_line = line[pos:].split()
                functional_info[try_str] = float(list_line[1])
        functional_str -= found
        if not functional_str: break
    return functional_info

# main program, get info from keyboard and convert the Landau parameters
if __name__ == '__main__':
    str_path = input('Please input path of thoout.dat:')
    input_type = input('Please input the type of input g (exp or not):')
    str_g = input('Please input the value of Landau parameter, use space to separate more than 1 input:')
    input_g = np.array([float(s) for s in str_g.split()])
    file_path = Path(str_path)
    if not str_path.endswith('/thoout.dat'):
        file_path = file_path / 'thoout.dat'
    with open(file_path, 'r') as f:
        functional_info = read_hfb(f)
    rho_nm = functional_info['RHO_NM=']
    smass_nm = functional_info['SMASS_NM=']
    hbzero = functional_info['hbzero=']
    if input_type.upper()=='EXP':
        print(landau_from_exp(input_g, rho_nm, smass_nm, hbzero))
    else:
        print(landau_to_exp(input_g, rho_nm, smass_nm, hbzero))
        
