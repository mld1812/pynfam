r"""
    Author: Tong Li, Mich State Univ, 2020-
"""

from pynfam import pynfam_mpi_calc
# from pynfam.utilities.mpi_utils import do_mpi
from pynfam.utilities import mpi_utils as mu
import os
import sys
import numpy as np
import pandas as pd
import functools
printf = functools.partial(print, flush=True) # force output flush
def printf_err(*msg): # function to print error message to both stdout and stderr
    printf(*msg, file=sys.stdout)
    printf(*msg, file=sys.stderr)
try:
    from pynfam_default import default_setts
except ImportError:
    printf('Warning: Cannot import default_setts from pynfam_default.py! Caution when using pynfam_fit_wrapper or pynfam_residual.')
    default_setts = {}
import pickle
from pynfam.config import DMASS_NH, DEFAULTS, EPSILON
import shutil
from copy import deepcopy
try:
    import requests
    import_requests = True
except ImportError:
    printf('Warning: Module requests cannot be imported. Assume mass file exists for GT/SD calcs.')
    import_requests = False
import io
import tarfile
import json
from scipy.optimize import minimize, least_squares
import datetime
import fnmatch

assert(sys.version_info[0] >= 3 and sys.version_info[1] >= 5) # Python >= 3.5 is required

# initialize a dict for transformation between one-level and two-level dicts
subkey_cat = {}
for key, subdict in DEFAULTS.items():
    for subdict2 in subdict.values():
        for final_key in subdict2.keys():
            subkey_cat[final_key] = key

# Landau parameters to use
subkey_cat['x1'] = 'scaled_params'
subkey_cat['g0p'] = 'scaled_params'
subkey_cat['g1p'] = 'scaled_params'
subkey_cat['h0p'] = 'scaled_params'
# scaled isoscalar pairing
subkey_cat['vpair_t0_scaled'] = 'scaled_params'
# pynfam_inputs
subkey_cat['gs_def_scan'] = 'pynfam_inputs'
subkey_cat['beta_type'] = 'pynfam_inputs'
subkey_cat['fam_ops'] = 'pynfam_inputs'


class energy_convert():
    # class for energy conversion between E_QRPA and E_x (excitation energy w.r.t. daughter g.s.)
    # see the appendix of PRC 101, 044305 (2020)

    def __init__(self, data_df):
        # read atomic mass info into a DataFrame
        mass_file = 'massround.mas20.txt'
        if (not os.path.isfile(mass_file)) and import_requests:
            url = 'https://www-nds.iaea.org/amdc/ame2020/massround.mas20.txt'
            headers = {'User-Agent': 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:77.0) Gecko/20100101 Firefox/77.0'}
            s = requests.get(url, headers = headers)
            with open(mass_file, 'wb') as file_object:
                file_object.write(s.content)
            mass_file = io.StringIO(s.text)
        self.df = pd.read_fwf(mass_file, skiprows=34, widths=[1,8,5,5,4,4,14,11,11,9,3,11,9,4,13,11], index_col=False, \
                              usecols=[1,2,4,6], names=['N','Z','element','mass_excess'], \
                              converters={'mass_excess': lambda x: np.nan if x.endswith('#') else float(x)})
        # drop unnecessary masses
        N_filter = data_df['N'].tolist() + (data_df['N']+1).tolist() + (data_df['N']-1).tolist()
        Z_filter = data_df['Z'].tolist() + (data_df['Z']+1).tolist() + (data_df['Z']-1).tolist()
        self.df.drop(self.df[(~self.df['N'].isin(N_filter)) | (~self.df['Z'].isin(Z_filter))].index, inplace=True)
        self.df['mass_excess'] /= 1000 # keV to MeV

    def _get_masses_(self, Z, N, beta_type):
        Delta_Tz = -1 if beta_type == '-' else +1
        m_p = self.df.loc[(self.df['Z'] == Z) & (self.df['N'] == N)]['mass_excess'].iat[0] # parent nucleus mass
        m_d = self.df.loc[(self.df['Z'] == Z-Delta_Tz) & (self.df['N'] == N+Delta_Tz)]['mass_excess'].iat[0] # daughter nucleus mass
        return Delta_Tz, m_p, m_d

    def shift_QRPA_to_Ex(self, E_QRPA, Z, N, lambda_n, lambda_p, beta_type):
        factor, m_p, m_d = self._get_masses_(Z, N, beta_type)
        return E_QRPA - m_d + m_p + factor * (lambda_n - lambda_p + DMASS_NH)
    
    def shift_Ex_to_QRPA(self, Ex, Z, N, lambda_n, lambda_p, beta_type):
        factor, m_p, m_d = self._get_masses_(Z, N, beta_type)
        return Ex + m_d - m_p - factor * (lambda_n - lambda_p + DMASS_NH)


def transform_upper(col):
    # helper function to transform all strings to upper case, except nucleus labels
    if col.name != 'label':
        return [s.upper() if isinstance(s, str) else s for s in col]
    else:
        return col

def check_pynfam_inputs(input1, input2, checktype):
    # helper function to check consistency of "pynfam_inputs"
    # to avoid possible conflicts with new rerun modes added to pynfam, we use large ints below
    if checktype == 'hfb':
        ignore_key = 'ignore_nonconv' # ignore the change of ignore_nonconv
        copy_input1 = {k: v for k, v in input1['hfb_mode'].items() if k != ignore_key}
        copy_input2 = {k: v for k, v in input2['hfb_mode'].items() if k != ignore_key}
        if copy_input1 != copy_input2:
            return 30 # hfb inconsistency
        else:
            return 0
    else:
        if input1['fam_mode'] != input2['fam_mode']:
            return 20 # fam inconsistency
        else:
            return 0

def check_override_settings(input1, input2, checktype):
    # helper function to check consistency of "override_settings"
    # to avoid possible conflicts with new rerun modes added to pynfam, we use large ints below
    if checktype == 'hfb':
        copy_input1 = input1['hfb']
        copy_input2 = input2['hfb']
        if copy_input1 != copy_input2:
            # special case: number of shells is changed
            special_key = 'number_of_shells'
            special_flag = (special_key in copy_input1.keys()) + (special_key in copy_input2.keys())
            if special_flag == 2:
                # both have special_key, compare corresponding values
                if copy_input1[special_key] != copy_input2[special_key]:
                    return 40 # hfb inconsistency, no restart from binaries
            elif special_flag == 1:
                # one has special_key but the other does not
                return 40 # hfb inconsistency, no restart from binaries
            return 30 # hfb inconsistency
        else:
            return 0
    else:
        if input1['fam'] != input2['fam'] or input1['ctr'] != input2['ctr']:
            return 20 # fam inconsistency
        else:
            return 0


def get_rerun(checktype, dir_path, category, pynfam_inputs, override_settings):
    # helper function to set rerun mode
    rerun = 0 # default
    pynfam_inputs_old = {}
    override_settings_old = {}
    pkl_path = os.path.join(dir_path, 'pynfam_inputs_'+category+'.pkl')
    if os.path.isfile(pkl_path):
        # check consistencies to set rerun
        rerun_mapping = {0:0, 20:'FAM', 30:'HFB', 40:'HFB_NORESTART'}
        with open(pkl_path, 'rb') as pklfile:
            pynfam_inputs_old = pickle.load(pklfile)
            override_settings_old = pickle.load(pklfile)
        c1 = check_pynfam_inputs(pynfam_inputs, pynfam_inputs_old, checktype)
        c2 = check_override_settings(override_settings, override_settings_old, checktype)
        rerun = rerun_mapping[max(c1, c2)]
    else:
        if os.path.isdir(dir_path):
            dir_list = os.listdir(dir_path)
        else:
            dir_list = None
        if dir_list and any((s.isdigit() or s == 'meta') and (os.path.isdir(os.path.join(dir_path, s))) for s in dir_list): 
            # rerun='HFB_NORESTART' if pkl does not exist but other outputs are found
            rerun = 'HFB_NORESTART'
            printf("Warning: Not found pkl file but output directories exist! Rerun=HFB_NORESTART is imposed.")
    if pynfam_inputs_old:
        return rerun, pynfam_inputs_old, override_settings_old
    else:
        return rerun, None, None


def move_or_del_file_dir(src, dst=None):
    # function to move file / dir from src to dst; dst will be overwritten
    # if only src is provided, src will be removed
    if not os.path.exists(src):
        return
    if dst is None:
        dst = src
    elif os.path.realpath(src) == os.path.realpath(dst):
        return
    if os.path.exists(dst):
        if os.path.isfile(dst):
            os.remove(dst)
        elif os.path.isdir(dst):
            shutil.rmtree(dst)
        else:
            raise IOError
    if src != dst:
        shutil.move(src, dst)

def files_rename_or_remove(dir_path, extension, suffix, untar=False):
    # rename or remove files with certain extensions, recursively in dir_path

    if untar:
        for label in os.listdir(dir_path):
            # untar all the tar files
            label_path = os.path.join(dir_path, label)
            if os.path.isfile(label_path) and label.endswith('.tar') and (not label.endswith('_og.tar')):
                with tarfile.open(label_path) as ftar:
                    for member in ftar:
                        ftar.extract(member, path=dir_path, set_attrs=False)
                new_path = os.path.join(dir_path, label[:-4] + '_og.tar')
                move_or_del_file_dir(label_path, new_path)

    # rename or remove files
    for label in os.listdir(dir_path):
        label_path = os.path.join(dir_path, label)
        if os.path.isfile(label_path):
            if label.endswith(extension):
                if suffix: # rename file by adding suffix
                    new_path = os.path.join(dir_path, label + suffix)
                    move_or_del_file_dir(label_path, new_path)
                else: # remove file
                    os.remove(label_path)
        elif os.path.isdir(label_path):
            files_rename_or_remove(label_path, extension, suffix, untar)
        else:
            raise IOError


def savepkl_and_clean(rerun, dir_path, category, restart, pynfam_inputs, override_settings, check):
    # save input parameters to pkl / json file, and rename existing output files to avoid possible restart problems

    if not os.path.isdir(dir_path):
        os.makedirs(dir_path)
        
    if check:
        pklfile = open(os.path.join(dir_path, 'pynfam_inputs_'+category+'_check.pkl'), 'wb')
        jsonfile = open(os.path.join(dir_path, 'pynfam_inputs_'+category+'_check.json'), 'w')
    else:
        pkl_path = os.path.join(dir_path, 'pynfam_inputs_'+category+'.pkl')
        json_path = os.path.join(dir_path, 'pynfam_inputs_'+category+'.json')
        pklfile = open(pkl_path, 'wb')
        jsonfile = open(json_path, 'w')
        # rename old output files when rerun == FAM, HFB or HFB_NORESTART
        if rerun == 'FAM':
            search_dirs = {'fam_soln', 'beta_soln'}
        elif (rerun == 'HFB') or (rerun == 'HFB_NORESTART'):
            search_dirs = {'hfb_soln', 'fam_soln', 'beta_soln'}
        else:
            search_dirs = None
        if search_dirs:
            for label in os.listdir(dir_path):
                label_path = os.path.join(dir_path, label)
                if not os.path.isdir(label_path): continue
                if label == 'meta':
                    files_rename_or_remove(label_path, ('.out', '.ctr', '.dat', '.in', '.py'), '.replaced', untar=False)
                else:
                    if not label.isdigit(): continue
                    for subdir in os.listdir(label_path):
                        subdir_path = os.path.join(label_path, subdir)
                        if (not subdir in search_dirs) or (not os.path.isdir(subdir_path)): continue
                        # rename existing output files in subdir_path, based on whether restarting from binaries
                        if restart:
                            files_rename_or_remove(subdir_path, ('.out', '.ctr', '.dat', '.in', '.log'), '.replaced', untar=True)
                        else:
                            files_rename_or_remove(subdir_path, ('.out', '.ctr', '.dat', '.in', '.log', '.bin', '.hel', '.tar'), \
                                                   '.replaced', untar=False)
                
    # write new inputs to files
    pickle.dump(pynfam_inputs, pklfile)
    pickle.dump(override_settings, pklfile)
    pklfile.close()
    json.dump((pynfam_inputs, override_settings), jsonfile, indent=4)
    jsonfile.close()


def read_hfb(dir_path, data_num):
    # read functional info, Fermi energies and pairing gaps from HFB calcs

    soln_dir = 'hfb_soln'
    output_filename = 'thoout.dat'
    other_filenames = ('hfbtho_NAMELIST.dat', 'hfbtho_output.hel')
    functional_str = {'RHO_NM=', 'SMASS_NM=', 'hbzero=', 'sigma=', \
                    'CJ(1)=', 'CpV0(0)=', 'CpV0(1)=', 't0=', 't1=', 't2=', 't3=', 'te=', 'to=', \
                    'use tensor terms:'}
    bool_functional_str = {'use tensor terms:'}
    functional_info = {}
    result_str = {'lambda', 'delta(n,p)'}
    result_info = {'lambda':[], 'delta(n,p)':[]}

    num = 0
    nonConv_labels = []
    all_OK = True
    for label in sorted(os.listdir(dir_path)): # go through all the directory labels
        if not label.isdigit(): continue
        label_path = os.path.join(dir_path, label)
        if not os.path.isdir(label_path): continue
        num += 1
        if num > data_num: break
        output_path = os.path.join(label_path, soln_dir, output_filename)
        if not (os.path.isfile(output_path) and \
           all(os.path.isfile(os.path.join(label_path, soln_dir, f)) for f in other_filenames)):
            # at least one hfb file is missing
            all_OK = False
            result_info['lambda'].append((np.nan,np.nan))
            result_info['delta(n,p)'].append((np.nan, np.nan))
            # for filename in (output_filename, *other_filenames):
            #     # create empty files which will be skipped by pynfam
            #     open(os.path.join(label_path, soln_dir, filename), 'w').close()
            continue
        
        with open(output_path, 'r') as output_file:
            temp_functional_str = deepcopy(functional_str)
            temp_functional_info = {}
            temp_result_str = deepcopy(result_str)
            conv_flag = False
            for line in output_file:
                # find functional info
                if temp_functional_str:
                    found = set()
                    for try_str in temp_functional_str:
                        pos = line.find(try_str)
                        if pos >= 0:
                            found.add(try_str)
                            list_line = line[(pos+len(try_str)):].split()
                            if try_str in bool_functional_str:
                                temp = (list_line[0].upper() in {'T','TRUE'})
                            else:
                                temp = float(list_line[0].rstrip(';'))
                                if np.abs(temp) < 1e-10: temp = 0
                            temp_functional_info[try_str] = temp
                    temp_functional_str -= found
                # find Fermi energies and pairing gaps
                for try_str in temp_result_str:
                    pos = line.find(try_str)
                    if pos >= 0:
                        list_line = line[pos:].split()
                        result_info[try_str].append((float(list_line[3]), float(list_line[4])))
                        temp_result_str.remove(try_str)
                        break
                # check convergence
                if line.find('iteration converged') >= 0:
                    conv_flag = True
            # check completeness and consistency
            if temp_functional_str:
                all_OK = False
            elif functional_info:
                assert(functional_info == temp_functional_info)
            else:
                functional_info = temp_functional_info
            if temp_result_str:
                all_OK = False
                for try_str in temp_result_str:
                    result_info[try_str].append((np.nan, np.nan))
            # record nonconvergence
            if not conv_flag:
                nonConv_labels.append(label)

    if num < data_num: # the number of dirs is smaller than the number of data
        all_OK = False
    if nonConv_labels: # print nonconvergent points
        printf("Warning: HFB run at", dir_path, "has non-convergent point(s) with following label(s).")
        printf(*nonConv_labels)
    lambda_n, lambda_p = zip(*result_info['lambda']) if result_info['lambda'] else ((), ())
    pair_gap_n, pair_gap_p = zip(*result_info['delta(n,p)']) if result_info['delta(n,p)'] else ((), ())
    return all_OK, functional_info, lambda_n, lambda_p, pair_gap_n, pair_gap_p


def to_couplings(override_settings, functional_info):
    # Convert Landau parameters to time-odd couplings, and scaled isoscalar pairing to the unscaled one.
    # Follow the convention of HFODD (p15-16 in arXiv:0909.3626 [nucl-th])
    # For tensor terms see PRC 85, 014326 (2012)
    # When a Landau parameter is not given, the corresponding coupling will be determined by
    # override_* or follow the PNFAM code if override_* does not exist.

    def get_fam_sett(key):
        if key in override_settings['fam'].keys():
            return override_settings['fam'][key]
        elif key in DEFAULTS['fam']['INTERACTION'].keys():
            return DEFAULTS['fam']['INTERACTION'][key]
        else:
            return None

    if 'scaled_params' not in override_settings.keys(): return
    if not override_settings['scaled_params']: return
    params_dict = deepcopy(override_settings['scaled_params'])

    # functional info
    rho_nm = functional_info['RHO_NM=']; inverse_effmass = functional_info['SMASS_NM=']
    hbzero = functional_info['hbzero=']; sigma = functional_info['sigma=']
    t0 = functional_info['t0=']; t1 = functional_info['t1=']; t2 = functional_info['t2=']; t3 = functional_info['t3=']
    te = functional_info['te=']; to = functional_info['to=']
    CJ = functional_info['CJ(1)=']; use_tensor_terms = functional_info['use tensor terms:']
    if not use_tensor_terms:
        te = 0; to = 0
    
    # Deal with isoscalar pairing scaled w.r.t. the abosolute value of isovector pairing
    vpair_t0_scaled = params_dict.pop('vpair_t0_scaled', None)
    if vpair_t0_scaled is not None:
        vpair_t1 = get_fam_sett('vpair_t1')
        if vpair_t1 is None:
            vpair_t1 = np.abs(np.average((functional_info['CpV0(0)='], functional_info['CpV0(1)='])))
        override_settings['fam']['vpair_t0'] = vpair_t0_scaled * vpair_t1

    # Start converting from Landau parameters to couplings
    x1 = params_dict.pop('x1', None)
    g0p = params_dict.pop('g0p', None)
    g1p = params_dict.pop('g1p', None)
    h0p = params_dict.pop('h0p', None)
    gauge = get_fam_sett('require_gauge_invariance')
    kf = (3/2 * np.pi**2 * rho_nm)**(1/3)
    kf_square = kf*kf
    inverse_kf_square = 1/kf_square # 1/k_F^2
    inverse_n0 = np.pi**(4/3) * hbzero * inverse_effmass / ((3/2 * rho_nm)**(1/3)) # 1/N_0

    # h_0^\prime
    if h0p is not None:
        cf_mid = h0p * inverse_n0 # = cf*k_F^2 / 3
        cf_over_3 = cf_mid * inverse_kf_square # = cf / 3
        cf = 3 * cf_over_3 # cf = 3*h0p/N_0/k_F^2
        override_settings['fam']['override_cf'] = cf
    else:
        cf = get_fam_sett('override_cf')
        if cf is None:
            if gauge:
                cf = 0.0
            else:
                cf = -3.0/8.0  * (te - to)
        cf_over_3 = cf / 3
        cf_mid = cf_over_3 * kf_square # cf*k_F^2 / 3
    
    # g_1^\prime
    if g1p is not None:
        temp = g1p * inverse_n0 / 2
        ct_mid =  - temp - cf_mid # = ct*k_F^2
        ct = - temp * inverse_kf_square - cf_over_3 # ct = -g1p/2/N_0/k_F^2 - cf/3
        override_settings['fam']['override_ct'] = ct
    else:
        ct = get_fam_sett('override_ct')
        if ct is None:
            if gauge:
                ct = -CJ
            else:
                ct = -1.0/16.0 * (t1 - t2 - 2.0*te + 2.0*to)
        ct_mid = ct * kf_square # ct*k_F^2
    
    # x_1 and g_0^\prime
    cs0 = get_fam_sett('override_cs0')
    csr = get_fam_sett('override_csr')
    if g0p is not None:
        cs_nm = g0p * inverse_n0 / 2 - ct_mid - cf_mid # cs_nm = g0p/2/N0 - ct*k_F^2 - cf*k_F^2 / 3
        if x1 is not None:
            cs0 = cs_nm * x1 # x1 = cs0 / cs_nm
            csr = (cs_nm - cs0) / (rho_nm**sigma) # cs_nm = cs0 + csr * rho_nm^sigma
        else:
            if (cs0 is None) and (csr is not None):
                cs0 = cs_nm - csr * (rho_nm**sigma)
            else:
                # ignore the overridden value of csr when cs0 is not None
                cs0 = -1/8 * t0 if cs0 is None else cs0
                csr = (cs_nm - cs0) / (rho_nm**sigma)
    else:
        if x1 is not None:
            if cs0 is not None:
                # ignore the overridden value of csr when cs0 is not None
                csr = cs0 / (rho_nm**sigma) * (1 - x1) / x1
            elif csr is not None:
                cs0 = csr * (rho_nm**sigma) * x1 / (1 - x1)
            else:
                # when both are None, first determine cs_nm from Skyrme functional and then use x1
                cs_nm = -1/8 * t0 - 1/48 * t3 * (rho_nm**sigma)
                cs0 = cs_nm * x1 # x1 = cs0 / cs_nm
                csr = (cs_nm - cs0) / (rho_nm**sigma) # cs_nm = cs0 + csr * rho_nm^sigma
    override_settings['fam']['override_cs0'] = cs0
    override_settings['fam']['override_csr'] = csr


def thieleInterpolator(x, y):
    r"""
    Returns a funcion which evaluates a rational function interpolation
    built from x and y using thiele continued fractions.
    This method produces a rational function with the degree of the numerator
    being either the same or one more than the degree of the denominator.

    Args:
        x (ndarray): The x values for the fit.
        y (ndarray): The y values for the fit.

    Returns:
        ufunc
    """
    # Calculate the continued fraction coefficients
    p = [[yi]*(len(y)-i) for i, yi in enumerate(y)]
    for i in range(len(p)-1):
        p[i][1] = (x[i] - x[i+1])/(p[i][0] - p[i+1][0])
    for i in range(2, len(p)):
        for j in range(len(p)-i):
            p[j][i] = (x[j]-x[j+i])/(p[j][i-1]-p[j+1][i-1]) + p[j+1][i-2]
    p0 = p[0]

    # Evaluate the continued fraction
    def t(xin):
        a = 0
        for i in range(len(p0)-1, 1, -1):
            a = ((xin - x[i-1])/(p0[i] - p0[i-2] + a))
        return y[0] + ((xin-x[0])/(p0[1]+a))

    # Check the interpolation passes though the provided points
    for xx, yy in zip(x, y):
        diff   = abs(yy - t(xx))
        mean   = 0.5*(abs(yy + t(xx)))
        reldiff = (diff/(mean+EPSILON))
        if min(reldiff,diff) > 1e-5:
            printf(u"Warning: ratint does not pass through points precisely.")
            printf(u"Relative diff: {:}, Abs diff: {:} at x={:}.".format(reldiff, diff, xx))

    # Return the evaluating function
    return t

def strength_interp(x_in, y_in):
    r"""
    Returns a funcion which evaluates a rational function interpolation
    for the complex strength function $S(\omega)$. It is assumed that 
    $S(\omega^*)=S^*(\omega)$, i.e. all the poles are on the real axis 
    with real-valued residuals
    
    Since $S(\omega)$ has the form of $\sum_{n}\frac{S_n}{\omega-E_{n}}$,
    the degree of its numerator is one less than the degree of the denominator,
    so here we use `thieleInterpolator` to interpolate $1/S(\omega)$.

    Args:
        x (ndarray): Values of complex energies $\omega$.
        y (ndarray): Values of complex strengths $S(\omega)$.

    Returns:
        ufunc
    """
    assert(len(x_in)==len(y_in))
    # S(w*) = S*(w*)
    x = np.dstack((x_in, np.conj(x_in))).flatten()
    y = np.dstack((y_in, np.conj(y_in))).flatten()
    def t(omega):
        return 1.0/thieleInterpolator(x, 1.0/y)(omega)
    return t


def detaildf_process_GTSD(output_df, category, label, column_name, converter, Z, N, lambda_n, lambda_p, beta_type, 
                          detail_dir_path, detail_fnmatch, use_ratinterp):
    energy = output_df['Re(EQRPA)'].to_numpy()
    hw = output_df['Im(EQRPA)'].to_numpy()
    assert(np.all(hw == hw[0])); hw = hw[0] # all the half widths should be equal (CONSTR)
    strength = output_df[column_name].to_numpy()
    energy_shifted = converter.shift_QRPA_to_Ex(energy, Z, N, lambda_n, lambda_p, beta_type)
    output_df.insert(0, 'Ex', energy_shifted)
    peak_pos = energy_shifted[np.argmax(strength)] # peak position
    num_error = np.nan

    # read and convert DataFrame to obtain complex strengths
    strength_list = []
    for detail_file in sorted(fnmatch.filter(os.listdir(detail_dir_path), detail_fnmatch)):
        detail_path = os.path.join(detail_dir_path, detail_file)
        if not os.path.isfile(detail_path): continue
        detail_df = pd.read_csv(detail_path, delim_whitespace=True, header=0, comment=u'#', index_col=0)
        op_name = detail_file[:-4].replace('-','_')
        output_df['S('+op_name+')'] = detail_df['Re(Strength)'] + 1j*detail_df['Im(Strength)']
        output_df['Conv('+op_name+')'] = (detail_df['Conv'] == 'Yes')
        k = int(op_name[-1])
        strength_arr =  output_df['S('+op_name+')'].to_numpy()
        strength_list.append((k, strength_arr))

    # Use a rational interpolator for complex strengths to find a better peak.
    # If use_ratinterp>1, numerical error caused by FAM convergence epsilon will also be computed.
    if use_ratinterp:
        eps_pos_list = [None]
        if use_ratinterp>1: # To compute numerical error caused by convergence epsilon
            num_error = 0.0
            eps_pos_list.extend(range(len(output_df)))
        for eps_pos in eps_pos_list:
            for factor in (1, 1j):
                temp_strength_list = []
                if eps_pos is None:
                    # To find the peak in the original data
                    if factor == 1j: continue
                    temp_strength_list = strength_list
                else:
                    # To obtain derivatives using finite difference for numerical error calculation
                    for k, strength_arr in strength_list:
                        temp_strength_arr = strength_arr.copy()
                        temp_strength_arr[eps_pos] += factor*output_df.at[eps_pos, 'Si']
                        temp_strength_list.append((k, temp_strength_arr))
                interp_list = [(k, strength_interp(energy_shifted+1j*hw, strength_arr)) \
                               for k, strength_arr in strength_list]
                printf('Scipy.minimize for', category, label,'starts at', datetime.datetime.now())
                opt_result = minimize(lambda x: -np.sum([(1.0 if k==0 else 2.0)*np.imag(interp(x+1j*hw)) for k, interp in interp_list]), \
                                      peak_pos, method='Nelder-Mead', options={'xatol':1e-14, 'fatol':1e-14, 'maxfev':1000})
                printf('Scipy.minimize for', category, label,'finishes at', datetime.datetime.now())
                if eps_pos is None:
                    printf('Scipy minimize summary for', category, label, ':')
                    printf(opt_result)
                    peak_pos = opt_result.x[0]
                    if not opt_result.success:
                        printf('Warning: Scipy minimize fails for', category, label)
                        # print some extra information for futher check
                        printf(opt_result.final_simplex[0].tolist(), opt_result.final_simplex[1].tolist())
                        printf(opt_result.final_simplex[0][1]-opt_result.final_simplex[0][0], \
                               opt_result.final_simplex[1][1]-opt_result.final_simplex[1][0])
                    if not ((peak_pos>=energy_shifted.min()) and (peak_pos<=energy_shifted.max())):
                        # Warning when the peak is outside the search region
                        printf('Warning:', category, label, 'has a peak outside the search region! The result may be inaccurate.')
                else:
                    num_error += (opt_result.x[0]-peak_pos)**2
        num_error = np.sqrt(num_error)
    
    return peak_pos, num_error


def detaildf_process_HL(output_df, spherical_flag=False, GA=None):
    # process detail DataFrame for the half life IN PLACE
    
    output_df.drop(index=['Total-GT', 'Forbidden-K=0', 'Forbidden-K=1', 'Forbidden-K=2', \
                   'Forbidden-J=0', 'Forbidden-J=1', 'Forbidden-J=2'], inplace=True, errors='ignore')
    if GA is not None:
        # recalculate Rate and NumErr_Rate with new GA value
        output_df['Rate(s^-1)'] = output_df['Rate_*GA^0(s^-1)'] + GA * output_df['Rate_*GA^1(s^-1)'] \
                                 + GA**2 * output_df['Rate_*GA^2(s^-1)']
        output_df['NumErr_Rate(s^-1)'] = np.sqrt(output_df['NumErrRate^2_*GA^0(s^-2)'] \
                                        + GA**2 * output_df['NumErrRate^2_*GA^2(s^-2)'] \
                                        + GA**4 * output_df['NumErrRate^2_*GA^4(s^-2)'])

    # calculate total rate and half life, with negative rates discarded
    allowed_names = ('Allowed-Fermi','Allowed-GT_K=0','Allowed-GT_K=1')
    ff_names = ('Forbidden-(J,K)=(0,0)','Forbidden-(J,K)=(1,0)','Forbidden-(J,K)=(1,1)',\
                'Forbidden-(J,K)=(2,0)','Forbidden-(J,K)=(2,1)','Forbidden-(J,K)=(2,2)')
    rows = ['Total', 'Total-Allowed', 'Total-Forbidden'] # rows to update afterwards
    if GA is not None:
        rows.extend(allowed_names + ff_names)
    if spherical_flag:
        # deal with spherical case: copy data of K=0 in columns given by "cols" to those of K!=0
        # the factor of 2 or 4 comes from "Tk" in shapeFactor/prepStrengths
        cols = ['Rate(s^-1)', 'Rate_*GA^0(s^-1)', 'Rate_*GA^1(s^-1)', 'Rate_*GA^2(s^-1)', \
                'NumErr_Rate(s^-1)', 'NumErrRate^2_*GA^0(s^-2)', 'NumErrRate^2_*GA^2(s^-2)', 'NumErrRate^2_*GA^4(s^-2)']
        factors = (2,2,2,2,2,4,4,4)
        for name in allowed_names + ff_names:
            if name.endswith(('Fermi', '0', '0)')):
                continue
            rows.append(name)
            if name.endswith(')'):
                output_df.loc[name, cols] = factors * output_df.loc[name[:-2] + '0)', cols]
            else:
                output_df.loc[name, cols] = factors * output_df.loc[name[:-1] + '0', cols]
    output_df.at['Total-Allowed', 'Rate(s^-1)'] = allowed_rate = \
        output_df.loc[output_df.index.isin(allowed_names) & (output_df['Rate(s^-1)']>0)]['Rate(s^-1)'].sum()
    output_df.at['Total-Forbidden', 'Rate(s^-1)'] = ff_rate = \
        output_df.loc[output_df.index.isin(ff_names) & (output_df['Rate(s^-1)']>0)]['Rate(s^-1)'].sum()
    output_df.at['Total', 'Rate(s^-1)'] = total_rate = allowed_rate + ff_rate

    # update other parts of output_df
    with np.errstate(divide=u'ignore', invalid=u'ignore'): # ignore 0/0 warning
        output_df.loc[rows, 'Half-Life(s)'] = np.log(2) / output_df.loc[rows, 'Rate(s^-1)']
        output_df.loc[rows, 'Log10(Half-Life)'] = np.log10(output_df.loc[rows, 'Half-Life(s)'])
        output_df.loc[rows, 'NumErr_Half-Life(s)'] = output_df.loc[rows, 'NumErr_Rate(s^-1)'] \
            * np.abs(output_df.loc[rows, 'Half-Life(s)'] / output_df.loc[rows, 'Rate(s^-1)'])
        output_df.loc[rows, 'NumErr_Log10(Half-Life)'] = output_df.loc[rows, 'NumErr_Rate(s^-1)'] \
            / np.abs(output_df.loc[rows, 'Rate(s^-1)']) / np.log(10)
        output_df.at[u'%FF', 'Rate(s^-1)'] = ff_rate / total_rate * 100 # percentage of first-forbidden rate


def df_func_residual(category_df, sigma_func_warning=True, category='the whole table'):
    # apply functions in column "func" onto model_raw and exp_raw and calculate residuals IN PLACE
    # numerical errors in column "num_error_raw" are also propagated
    # results are stored in column "model", "exp" and "num_error"
    np.seterr('ignore')
    for col in ('sigma', 'func'):
        if col not in category_df: # fill nan if sigma / func does not exist
            category_df[col] = np.nan
            continue
        vals = category_df[col].to_numpy()
        if (not ((vals[0] == vals).all() or pd.isnull(vals).all())) and sigma_func_warning:
            # check if sigma / func are the same for one category
            printf('Warning:', category, col, 'values are not equal!')
    if 'exp_raw' not in category_df:
        if 'exp' in category_df: # save exp to exp_raw
            category_df['exp_raw'] = category_df['exp']
            printf('Warning:', category, 'has col "exp" moved to "exp_raw", causing inconsistency!')
        else:
            category_df['exp_raw'] = np.nan
    fun_model_results = []; fun_exp = []; fun_num_errors = []
    for ii, func in enumerate(category_df['func']):
        # apply func onto model and exp, and compute corresponding num_error
        fun = lambda x: x
        derivative  = lambda x: 1
        if isinstance(func, str):
            if func == 'LOG' or func == 'LOG10':
                fun = np.log10
                derivative = lambda x: 1.0/x/np.log(10)
            elif func == 'LN':
                fun = np.log
                derivative = lambda x: 1.0/x
            elif func == 'RATE':
                fun = lambda x: np.log(2)/x
                derivative = lambda x: -np.log(2)/x/x
            else:
                if sigma_func_warning:
                    printf('Warning: Unsupported func', func, 'found in category', category, '! Fall back to identity function.')
        elif pd.isnull(func):
            pass
        else:
            if sigma_func_warning:
                printf('Warning: Unsupported func found in category', category, '! Fall back to identity function.')
        fun_model_results.append(fun(category_df['model_raw'].iat[ii]))
        fun_exp.append(fun(category_df['exp_raw'].iat[ii]))
        fun_num_errors.append(np.abs(derivative(category_df['model_raw'].iat[ii]))*category_df['num_error_raw'].iat[ii])
    category_df['exp'] = fun_exp
    category_df['model'] = fun_model_results
    category_df['num_error'] = fun_num_errors
    category_df['residual'] = category_df['model'] - category_df['exp']
    category_df['norm_res'] = category_df['residual'] / category_df['sigma']


def df_clean(output_data, exp_raw_in_input=True):
    # drop columns with all NaN, plus unnecessary columns when "func" does not exist, IN PLACE
    output_data.dropna(axis=1, how='all', inplace=True) # remove columns with all NaN
    if 'func' not in output_data: # drop columns related to the application of func
        cols_to_drop = ['model_raw', 'num_error_raw']
        if not exp_raw_in_input:
            cols_to_drop.append('exp_raw')
        output_data.drop(columns=cols_to_drop, inplace=True, errors='ignore')


def scratch_from_backup(backup_subdir_path, scratch_dir_path, categories):
    # start from a given backup dir; return whether it succeeds
    if os.path.isdir(backup_subdir_path):
        try:
            for item in os.listdir(backup_subdir_path): # copy from backup
                if item not in categories: continue
                item_path = os.path.join(backup_subdir_path, item)
                if not os.path.isdir(item_path): continue
                temp_path = os.path.join(scratch_dir_path, item)
                # rename existing dirs/files as "_replaced"
                move_or_del_file_dir(temp_path, os.path.join(scratch_dir_path, item + '_replaced'))
                # copy; ignore the tar file which will be untarred later
                ftar_name = item + '.tar'
                shutil.copytree(item_path, temp_path, ignore=shutil.ignore_patterns(ftar_name))
                # untar files
                ftar_path = os.path.join(item_path, ftar_name)
                if os.path.isfile(ftar_path):
                    with tarfile.open(ftar_path) as ftar:
                        for member in ftar:
                            ftar.extract(member, path=temp_path, set_attrs=False)
            printf('Info: Start from backup', backup_subdir_path, 'succeeds.')
            flag = True
        except:
            printf('Warning: Start from backup', backup_subdir_path, 'fails!')
            flag = False
    else:
        printf('Info: Backup', backup_subdir_path, 'is not found. Start from current scratch.')
        flag = False
    return flag


def backup_num_limit(backup_max, backup_dir_path):
    # remove old backups to limit the total number of backups; return next backup subdir name
    backup_subdirs = sorted([int(x) for x in os.listdir(backup_dir_path) if x.isdigit()], reverse=True)
    if backup_subdirs:
        max_dir = backup_subdirs[0]
    else:
        max_dir = -1
    if backup_max > 0: # remove dirs to keep the number of existing backups below backup_max
        while len(backup_subdirs) >= backup_max:
            try:
                path_remove = os.path.join(backup_dir_path, str(backup_subdirs[-1]))
                move_or_del_file_dir(path_remove)
            except:
                printf("Warning: Backup", str(backup_subdirs[-1]), "cannot be removed! Skipped.")
                break
            backup_subdirs.pop()
    return str(max_dir+1)


def scratch_to_backup(backup_subdir_path, scratch_dir_path, ignore_tuple, categories, result_filename, overwrite_flag=False):
    # copy scratch to backup
    ignore_func = shutil.ignore_patterns(*ignore_tuple)
    try:
        # remove existing directory/file
        if os.path.exists(backup_subdir_path):
            move_or_del_file_dir(backup_subdir_path)
        # copy directories
        for item in os.listdir(scratch_dir_path):
            if item not in categories: continue
            dir_path = os.path.join(scratch_dir_path, item)
            if not os.path.isdir(dir_path): continue
            dst_dir_path = os.path.join(backup_subdir_path, item)
            # copy
            shutil.copytree(dir_path, dst_dir_path, ignore=ignore_func)
            # tar dirs in dst_dir_path to reduce file numbers
            with tarfile.open(os.path.join(dst_dir_path, item) +'.tar', 'w') as ftar:
                for dir_to_tar in os.listdir(dst_dir_path):
                    dir_to_tar_path = os.path.join(dst_dir_path, dir_to_tar)
                    if os.path.isdir(dir_to_tar_path):
                        ftar.add(dir_to_tar_path, arcname=os.path.relpath(dir_to_tar_path, dst_dir_path))
                        shutil.rmtree(dir_to_tar_path)
        # copy result files
        shutil.copy2(os.path.join(scratch_dir_path, result_filename +'.pkl'), backup_subdir_path)
        shutil.copy2(os.path.join(scratch_dir_path, result_filename +'.csv'), backup_subdir_path)
        if overwrite_flag:
            printf("Info: Backup to existing directory", backup_subdir_path, "succeeds.")
        else:
            printf("Info: Backup to directory", backup_subdir_path, "succeeds.")
    except:
        printf("Warning: Backup to directory", backup_subdir_path, "fails!")


def pynfam_fit_wrapper_nonroot(comm, root_rank):
    # Work function for non-root processes
    assert(mu.do_mpi)
    val = None
    while True:
        val = comm.bcast(val, root=root_rank) # receive parameters from root
        if isinstance(val, tuple):
            # run pynfam
            pynfam_mpi_calc(*val, comm)
        else:
            # return what's broadcasted
            return val
            

def pynfam_fit_wrapper_root(input_data, categories, override_setts_fit, override_setts_cat, return_df, clean_files,
                            exes_dir, scratch_dir, backup_dir, backup_option, backup_max, start_from_backup,
                            nr_parallel_calcs, nthreads_per_calc, rerun_mode, ignore_hfb_nonconv, hfb_only, use_ratinterp,
                            comm, rank, check, result_filename, **kwargs):
    # Work function for root process

    # check / read input_data
    if isinstance(input_data, str):
        input_data = pd.read_csv(input_data, header=0, comment='#')
    else:
        assert(isinstance(input_data, pd.DataFrame))
    
    # tidy up input_data and sort it by categories according to the order of their first appearances
    all_data = input_data.apply(transform_upper) # change to upper case
    all_categories = tuple(all_data['category'])
    all_categories_pure = tuple(sorted(set(all_categories), key=all_categories.index)) # remove duplicates but keep order
    all_data['cat_index'] = tuple(map(all_categories_pure.index, all_categories)) # add auxialiary column for sort
    all_data.sort_values(by='cat_index', inplace=True, kind='mergesort', ascending=True)
    all_data.drop('cat_index', axis=1, inplace=True) # drop auxiliary column

    # check list / tuple inputs, or convert dict inputs to tuples and check
    if categories is None:
        categories = all_categories_pure
    if not isinstance(categories, (list, tuple)):
        categories = (categories,)
    categories = tuple(map(str.upper, categories)) # change to upper case
    def dict_to_tuple_by_cat(d, cats): # convert dict to tuple according to categories
        if not isinstance(d, dict):
            return (d,) * len(categories)
        temp_d = {str(k).upper(): v for k, v in d.items()}
        return tuple(temp_d.get(cat, None) for cat in cats)
    nr_parallel_calcs = dict_to_tuple_by_cat(nr_parallel_calcs, categories)
    nthreads_per_calc = dict_to_tuple_by_cat(nthreads_per_calc, categories)
    rerun_mode = tuple(x.upper() if isinstance(x, str) else x \
                       for x in dict_to_tuple_by_cat(rerun_mode, categories)) # change to upper case
    assert(all(x in {None,0,1,'FAM','HFB','HFB_NORESTART'} for x in rerun_mode))
    ignore_hfb_nonconv = tuple(x if x is not None else 2 for x in \
                            dict_to_tuple_by_cat(ignore_hfb_nonconv, categories)) # change None to 2
    hfb_only = dict_to_tuple_by_cat(hfb_only, categories)

    # check other inputs
    assert(isinstance(override_setts_fit, dict))
    assert(isinstance(override_setts_cat, dict))
    assert(use_ratinterp in {0,1,2})
    assert(backup_option in {0,1,2,-1,-2})
    assert(isinstance(backup_max, int))
    assert(isinstance(backup_dir, str))
    assert(isinstance(start_from_backup, int))
    assert(isinstance(scratch_dir, str))
    assert(isinstance(exes_dir, str))

    # prepare scratch dir
    scratch_dir_path = os.path.realpath(scratch_dir)
    if not os.path.isdir(scratch_dir_path):
        os.makedirs(scratch_dir_path)

    # prepare backup dir and get backup files if we need
    if backup_option or start_from_backup >= 0:
        backup_dir_path = os.path.realpath(backup_dir)
        if not os.path.isdir(backup_dir_path):
            os.makedirs(backup_dir_path)
        if start_from_backup >= 0:
            if not scratch_from_backup(os.path.join(backup_dir_path, str(start_from_backup)), scratch_dir_path, categories):
                start_from_backup = -1
    
    # main work starts
    support_categories = {'GT', 'SD', 'HL', 'GAP'}
    category_set = set(categories) & set(all_categories_pure) & support_categories
    all_model_results = []
    all_rerun_equal_zero = True
    if ('GT' in category_set) or ('SD' in category_set): # prepare for energy conversion
        converter = energy_convert(all_data.loc[(all_data['category']=='GT') | (all_data['category']=='SD')])
    # rename existing result files
    move_or_del_file_dir(os.path.join(scratch_dir_path, result_filename +'.pkl'), \
                         os.path.join(scratch_dir_path, result_filename +'.pkl.replaced'))
    move_or_del_file_dir(os.path.join(scratch_dir_path, result_filename +'_other.pkl'), \
                         os.path.join(scratch_dir_path, result_filename +'_other.pkl.replaced'))
    move_or_del_file_dir(os.path.join(scratch_dir_path, result_filename +'.csv'), \
                         os.path.join(scratch_dir_path, result_filename +'.csv.replaced'))

    for ind, category in enumerate(categories):
        printf('Category:', category)
        model_results = []
        output_dfs = []
        conv_info = []
        maxsi_info = []
        num_errors = []

        # avoid duplicate / unsupported category or category with no data
        if category not in category_set:
            printf('Warning: Category', category, 'is unsupported, missing in input data or duplicate! Skipped.')
            continue
        category_set.remove(category)
        data = all_data.loc[all_data['category']==category] # get data for a given category

        # construct pynfam_inputs and override_settings for pynfam, hfb only
        # start from pynfam_inputs, whose default values are hard coded here
        pynfam_inputs = {
        'directories': {
            'outputs' : category,
            'exes'    : exes_dir,
            'scratch' : scratch_dir
            },
        'nr_parallel_calcs': nr_parallel_calcs[ind],
        'nthreads_per_calc': nthreads_per_calc[ind],
        'hfb_mode': {
            'dripline_mode'  : 0,
            'ignore_nonconv' : ignore_hfb_nonconv[ind],
            'gs_def_scan'    : (-2, (-0.2,0.0,0.2)),
            'retry_nonconv'  : True,
            'extra_inputs'   : ()
            },
        'fam_mode': {
            'fam_contour'       : None,
            'beta_type'         : '-',
            'fam_ops'           : 'ALL'
            }
        }
        # override_settings from various dicts
        override_settings = {'pynfam_inputs':{}, 'hfb': {}, 'fam': {}, 'ctr': {}, 'psi': {}, 'scaled_params': {}}
        for param_dict in (default_setts, override_setts_cat.get(category,{}), override_setts_fit, kwargs):
            if not (param_dict.keys() & override_settings.keys()):
                # one-level dict input, transform to two-level
                temp_dict = {'hfb':{}, 'fam':{}, 'ctr': {}, 'psi': {}, 'scaled_params': {}}
                for sub_key, val in param_dict.items():
                    if sub_key in subkey_cat.keys():
                        temp_dict[subkey_cat[sub_key]][sub_key] = val
                    else:
                        printf('Warning: One-level setting', sub_key, 'is ignored!')
            else:
                # two-level dict input
                temp_dict = param_dict
            for key, sub_dict in temp_dict.items():
                if key not in override_settings.keys():
                    printf('Warning: First-level setting', key, 'is ignored!')
                    continue
                for sub_key, val in sub_dict.items():
                    if sub_key in subkey_cat.keys():
                        override_settings[key][sub_key] = val
                    else:
                        printf('Warning: Second-level setting', sub_key, 'is ignored!')
        # update pynfam_inputs from override_settings['pynfam_inputs']
        pynfam_inputs_update = override_settings.pop('pynfam_inputs')
        if 'gs_def_scan' in pynfam_inputs_update:
            pynfam_inputs['hfb_mode']['gs_def_scan'] = pynfam_inputs_update['gs_def_scan']
        if 'beta_type' in pynfam_inputs_update:
            pynfam_inputs['fam_mode']['beta_type'] = pynfam_inputs_update['beta_type']
        if 'fam_ops' in pynfam_inputs_update:
            pynfam_inputs['fam_mode']['fam_ops'] = pynfam_inputs_update['fam_ops']
        # neutron / proton numbers
        override_settings['hfb']['proton_number'] = tuple(data['Z'])
        override_settings['hfb']['neutron_number'] = tuple(data['N'])
        # deformation info given in input dataframe
        if 'deformation' in data.columns:
            hfb_gs_def_scan = []
            fail_row = []
            for row, deformation in zip(data.index, data['deformation']):
                if isinstance(deformation, str):
                    scan_val = []
                    if deformation.find('O') >= 0: # oblate
                        scan_val.append(-0.2)
                    if deformation.find('S') >= 0: # spherical
                        scan_val.append(0.0)
                    if deformation.find('P') >= 0: # prolate
                        scan_val.append(0.2)
                    if not scan_val: # deal with string input like '(0.0,0.2)', '[0.0,0.2]', '0.0,0.2'
                        try:
                            scan_val = map(float, deformation.strip('()[] ').split(','))
                        except:
                            fail_row.append(row)
                            scan_val = (-0.2,0.0,0.2)
                    hfb_gs_def_scan.append((-2, tuple(scan_val)))
                elif isinstance(deformation, tuple):
                    hfb_gs_def_scan.append(deformation)
                else:
                    fail_row.append(row)
                    hfb_gs_def_scan.append(pynfam_inputs['hfb_mode']['gs_def_scan']) # default value
            if fail_row:
                printf('Info: Deformation scan uses default', pynfam_inputs['hfb_mode']['gs_def_scan'], 'in row(s)', fail_row)
            if len(set(hfb_gs_def_scan)) == 1:
                pynfam_inputs['hfb_mode']['gs_def_scan'] = hfb_gs_def_scan[0]
            else:
                pynfam_inputs['hfb_mode']['gs_def_scan'] = hfb_gs_def_scan
        else:
            hfb_gs_def_scan = [pynfam_inputs['hfb_mode']['gs_def_scan']] * len(data)
            printf('Info: Deformation scan uses default', pynfam_inputs['hfb_mode']['gs_def_scan'])

        # check input and set rerun mode
        dir_path = os.path.join(scratch_dir_path, category)
        if rerun_mode[ind] is None:
            rerun, pynfam_inputs_hfb, override_settings_hfb = get_rerun('hfb', dir_path, category, pynfam_inputs, override_settings)
            if pynfam_inputs_hfb:
                # update pynfam_inputs_hfb and override_settings_hfb, except fam/ctr part; both dicts will be saved in the pkl file
                for key in pynfam_inputs.keys():
                    if key != 'fam_mode':
                        pynfam_inputs_hfb[key] = pynfam_inputs[key]
                for key in override_settings.keys():
                    if key != 'fam' and key != 'ctr':
                        override_settings_hfb[key] = override_settings[key]
            else:
                pynfam_inputs_hfb, override_settings_hfb = pynfam_inputs, override_settings
        else:
            rerun = rerun_mode[ind]
            pynfam_inputs_hfb, override_settings_hfb = pynfam_inputs, override_settings
            if rerun in {0,1,'FAM'}:
                printf("Warning: It's dangerous to set rerun_mode=0, 1 or FAM without checking the pkl file for the HFB part.")
        if 'restart_file' not in override_settings['hfb']:
            override_settings['hfb']['restart_file'] = DEFAULTS['hfb']['HFBTHO_ITERATIONS']['restart_file']
        if rerun == 'HFB_NORESTART': # when number of shells is changed, hfb cannot restart from existing binary files.
            override_settings['hfb']['restart_file'] = override_settings_hfb['hfb']['restart_file'] = 1
        pynfam_inputs['rerun_mode'] = pynfam_inputs_hfb['rerun_mode'] = rerun if isinstance(rerun, int) else 0
        printf('HFB rerun_mode is', rerun)
        all_rerun_equal_zero = all_rerun_equal_zero and (rerun == 0)

        # save pkl and file clean before starting the main calculation
        savepkl_and_clean(rerun, dir_path, category, override_settings['hfb']['restart_file']<0, \
                          pynfam_inputs_hfb, override_settings_hfb, check)
        # rename existing result files
        move_or_del_file_dir(os.path.join(dir_path, result_filename + '_' + category + '.pkl'), \
                             os.path.join(dir_path, result_filename + '_' + category + '.pkl.replaced'))
        move_or_del_file_dir(os.path.join(dir_path, result_filename + '_' + category + '.csv'), \
                             os.path.join(dir_path, result_filename + '_' + category + '.csv.replaced'))

        # read info from hfb calcs to check whether we can avoid calling pynfam_mpi_calc
        all_hfb_OK = False
        if rerun in {0,1,'FAM'}:
            all_hfb_OK, functional_info, lambda_n, lambda_p, pair_gap_n, pair_gap_p = read_hfb(dir_path, len(data))

        # run pynfam, hfb only
        if rerun != 1 and (not all_hfb_OK):
            # when rerun==1 or hfb info is extracted successfully, no HFB calc will be performed
            if mu.do_mpi:
                pynfam_inputs, override_settings, check = comm.bcast((pynfam_inputs, override_settings, check), root=rank)
            printf('Call of pynfam_mpi_calc (HFB part) for', category, 'starts at', datetime.datetime.now())
            pynfam_mpi_calc(pynfam_inputs, override_settings, check, comm)
            printf('Call of pynfam_mpi_calc (HFB part) for', category, 'finishes at', datetime.datetime.now())
            # read info from hfb calcs (again)
            all_hfb_OK, functional_info, lambda_n, lambda_p, pair_gap_n, pair_gap_p = read_hfb(dir_path, len(data))
        else:
            printf("Call of pynfam_mpi_calc (HFB part) is skipped for", category)

        # Shall we skip FAM part?
        if not all_hfb_OK:
            if check:
                printf("Warning: PyNFAM run of", category, "does not have valid HFB files to check FAM part! FAM part is skipped.")
            else:
                printf_err("Error: PyNFAM run of", category, "does not give valid HFB files! FAM part is skipped.")
        elif hfb_only[ind]:
            printf("Call of pynfam_mpi_calc (FAM part) is skipped for", category, 'because of hfb_only setting.')
        elif category == 'GAP':
            printf("Call of pynfam_mpi_calc (FAM part) is skipped for", category)
        if (not all_hfb_OK) or (hfb_only[ind] and category != 'GAP'):
            # when category is 'GAP', we still need to run codes below to collect results
            all_model_results.append(deepcopy(data))
            continue

        # beta type
        if 'beta_type' in data.columns:
            btypes = [b if not pd.isnull(b) else pynfam_inputs['fam_mode']['beta_type'] for b in data['beta_type']]
            if len(set(btypes)) == 1:
                pynfam_inputs['fam_mode']['beta_type'] = btypes[0]
            else:
                pynfam_inputs['fam_mode']['beta_type'] = btypes
        else:
            btypes = (pynfam_inputs['fam_mode']['beta_type'],) * len(data)
        
        # construct input for pynfam, fam only
        if category == 'GT' or category == 'SD': # Gamow-Teller or spin-dipole resonance
            if category == 'GT':
                long_category = 'GAMOWTELLER'; kval_category = 'GT'
            else:
                long_category = kval_category = 'SPINDIPOLE'
            to_couplings(override_settings, functional_info) # convert scaled parameters to couplingss
            pynfam_inputs['fam_mode']['fam_contour'] = 'CONSTR'
            # energy conversion, see the appendix of PRC 101, 044305 (2020)
            override_settings['ctr']['energy_min'] = tuple(converter.shift_Ex_to_QRPA(Emin, Z, N, ln, lp, btype) \
                                                           for Emin, Z, N, ln, lp, btype in \
                                                           zip(data['energy_min'], data['Z'], data['N'], lambda_n, lambda_p, btypes))
            override_settings['ctr']['energy_max'] = tuple(converter.shift_Ex_to_QRPA(Emax, Z, N, ln, lp, btype) \
                                                           for Emax, Z, N, ln, lp, btype in \
                                                           zip(data['energy_max'], data['Z'], data['N'], lambda_n, lambda_p, btypes))
            if 'half_width' in data.columns:
                if 'half_width' not in override_settings['ctr']:
                    override_settings['ctr']['half_width'] = DEFAULTS['CONSTR']['half_width']
                hw_tuple = tuple(x if not np.isnan(x) else override_settings['ctr']['half_width'] for x in data['half_width'])
                override_settings['ctr']['half_width'] = hw_tuple
            if 'de_hw_ratio' in data.columns:
                if 'de_hw_ratio' not in override_settings['ctr']:
                    override_settings['ctr']['de_hw_ratio'] = DEFAULTS['CONSTR']['de_hw_ratio']
                rt_tuple = tuple(x if not np.isnan(x) else override_settings['ctr']['de_hw_ratio'] for x in data['de_hw_ratio'])
                override_settings['ctr']['de_hw_ratio'] = rt_tuple
            if 'deformation' in data.columns: # reduce computation for spherical case
                ops = []
                sph_row = []
                for row, deformation in zip(data.index, hfb_gs_def_scan):
                    if deformation == (-2, (0.0,)): # spherical case
                        ops.append((kval_category, 0))
                        sph_row.append(row)
                    else:
                        ops.append(long_category)
                if sph_row:
                    printf('Info:', category, 'operator uses K=0 only for spherical system(s) in row(s)', sph_row)
                pynfam_inputs['fam_mode']['fam_ops'] = ops[0] if len(set(ops)) == 1 else ops
            else:
                pynfam_inputs['fam_mode']['fam_ops'] = long_category
            # remove unnecessary ctr keys
            override_settings['ctr'] = {x: val for x, val in override_settings['ctr'].items() \
                                       if x in DEFAULTS['ctr']['CONSTR'].keys() or x == 'energy_min' or x == 'energy_max'}
        elif category == 'HL': # beta decay half-life
            to_couplings(override_settings, functional_info) # convert scaled parameters to couplings
            pynfam_inputs['fam_mode']['fam_contour'] = 'CIRCLE'
            if 'op' in data.columns: # fam operators
                ops_str = []
                fail_row = []
                for row, operator in zip(data.index, data['op']):
                    if operator in {'ALL','ALLOWED','FORBIDDEN','GAMOWTELLER','SPINDIPOLE','0+','1+','0-','1-','2-'}:
                        ops_str.append(operator)
                    else:
                        ops_str.append(pynfam_inputs['fam_mode']['fam_ops']) # default value
                        fail_row.append(row)
                if fail_row:
                    printf('Info: HL operator uses default', pynfam_inputs['fam_mode']['fam_ops'], 'in row(s)', fail_row)
            else:
                ops_str = [pynfam_inputs['fam_mode']['fam_ops']] * len(data) # default value
                printf('Info: HL operator uses default', pynfam_inputs['fam_mode']['fam_ops'])
            if 'deformation' in data.columns: # reduce computation for spherical case
                ops = []
                sph_row = []
                for row, deformation, op_str in zip(data.index, hfb_gs_def_scan, ops_str):
                    if deformation == (-2, (0.0,)): # spherical case
                        ops.append((op_str, 0))
                        sph_row.append(row)
                    else:
                        ops.append(op_str)
                if sph_row:
                    printf('Info: HL operator uses K=0 only for spherical system(s) in row(s)', sph_row)
            else:
                ops = ops_str
            pynfam_inputs['fam_mode']['fam_ops'] = ops[0] if len(set(ops)) == 1 else ops
            # remove unnecessary ctr keys
            override_settings['ctr'] = {x: val for x, val in override_settings['ctr'].items() if x in DEFAULTS['ctr']['CIRCLE'].keys()}
        elif category == 'GAP': # pairing gap (HFB only, no FAM to be done)
            pass
        else:
            continue # unsupported category, ignore

        if pynfam_inputs['fam_mode']['fam_contour'] is not None:
            # check input and set rerun mode (again)
            if rerun_mode[ind] is None:
                if (rerun !='HFB') and (rerun != 'HFB_NORESTART'):
                    # if hfb is not rerun, check fam input; otherwise, fam will be rerun without checking
                    rerun, _, _ = get_rerun('fam', dir_path, category, pynfam_inputs, override_settings)
            else:
                rerun = rerun_mode[ind]
                if rerun == 0 or rerun == 1:
                    printf("Warning: It's dangerous to set rerun_mode=0 or 1 without checking the pkl file for the FAM part.")
            if 'use_fam_storage' not in override_settings['fam']:
                override_settings['fam']['use_fam_storage'] = DEFAULTS['fam']['GENERAL']['use_fam_storage']
            if (rerun == 'HFB_NORESTART') and (override_settings['fam']['use_fam_storage'] == -1):
                # when number of shells is changed, fam cannot restart from existing binary files.
                override_settings['fam']['use_fam_storage'] = 1
            pynfam_inputs['rerun_mode'] = rerun if isinstance(rerun, int) else 0
            printf('FAM rerun_mode is', rerun)
            all_rerun_equal_zero = all_rerun_equal_zero and (rerun == 0)

            # save pkl and file clean before start the main calculation (again)
            if (rerun == 'HFB') or (rerun == 'HFB_NORESTART'): rerun = 'FAM'
            savepkl_and_clean(rerun, dir_path, category, override_settings['fam']['use_fam_storage']<0, \
                              pynfam_inputs, override_settings, check)

            # run pynfam, fam only
            if mu.do_mpi:
                pynfam_inputs, override_settings, check = comm.bcast((pynfam_inputs, override_settings, check), root=rank)
            printf('Call of pynfam_mpi_calc (FAM part) for', category, 'starts at', datetime.datetime.now())
            pynfam_mpi_calc(pynfam_inputs, override_settings, check, comm)
            printf('Call of pynfam_mpi_calc (FAM part) for', category, 'finishes at', datetime.datetime.now())
        
        if check: continue # finish check

        # get output from files
        soln_dir = 'beta_soln'
        detail_dir = 'fam_soln'
        if category == 'GT': # Gamow-Teller resonance
            output_filename = 'totalGT_str.out'
            column_name = 'Total-GT'
            detail_fnmatch = 'GT-K[0-1].out'
        elif category == 'SD': # spin-dipole resonance
            output_filename = 'totalSD_str.out'
            column_name = 'Total-SD'
            detail_fnmatch = 'RS[0-2]-K[0-2].out'
        elif category == 'HL': # beta decay half-life
            output_filename = 'beta.out'
            column_name = 'Half-Life(s)'
        elif category == 'GAP': # pairing gap
            soln_dir = 'hfb_soln'
            detail_dir = 'hfb_soln'
            output_filename = 'thoout.dat'
        else:
            continue # unsupported category, ignore
        
        num = -1
        for label in sorted(os.listdir(dir_path)): # go through all the directory labels
            if not label.isdigit(): continue
            label_path = os.path.join(dir_path, label)
            if not os.path.isdir(label_path): continue
            num += 1
            if num >= len(data): break
            output_path = os.path.join(label_path, soln_dir, output_filename)
            detail_dir_path = os.path.join(label_path, detail_dir)
            if (not os.path.isfile(output_path)) or (not os.path.isdir(detail_dir_path)):
                printf_err("Error: PyNFAM run of", category, label, "does not give necessary output file(s)!")
                model_results.append(np.nan); output_dfs.append(np.nan); num_errors.append(np.nan)
                continue

            if category == 'GT' or category == 'SD': # Gamow-Teller or spin-dipole resonance
                output_df = pd.read_csv(output_path, delim_whitespace=True, header=0, comment=u'#', index_col=0)
                peak_pos, num_error = detaildf_process_GTSD(output_df, category, label, column_name, converter, \
                                        data['Z'].iat[num], data['N'].iat[num], lambda_n[num], lambda_p[num], btypes[num], \
                                        detail_dir_path, detail_fnmatch, use_ratinterp)
                model_results.append(peak_pos)
                num_errors.append(num_error)

            elif category == 'HL': # beta decay half-life
                output_df = pd.read_csv(output_path, delim_whitespace=True, header=0, comment=u'#', index_col=0)
                detaildf_process_HL(output_df, hfb_gs_def_scan[num] == (-2, (0.0,)))
                model_results.append(output_df.at['Total', column_name])
                num_errors.append(output_df.at['Total', 'NumErr_'+column_name])

            elif category == 'GAP': # pairing gap
                if 'op' in data.columns:
                    select = data['op'].iat[num]
                else:
                    select = None
                if select == 'N' or select == 'NEUTRON':
                    model_results.append(pair_gap_n[num])
                elif select == 'P' or select == 'PROTON':
                    model_results.append(pair_gap_p[num])
                else:
                    model_results.append((pair_gap_n[num], pair_gap_p[num]))
                output_df = np.nan
                num_errors.append(np.nan)

            output_dfs.append(output_df)
        # loop over labels ends here

        # read convergence info and max si
        log_path = os.path.join(dir_path, 'meta', 'logfile_master.dat')
        if not os.path.isfile(log_path):
            printf_err("Error: PyNFAM run of", category, "does not give a log file! Nonconvergence assumed.")
            conv_info.extend((False,)*len(data))
            maxsi_info.extend((np.nan,)*len(data))
        else:
            log_df = pd.read_csv(log_path, delim_whitespace=True, header=0, comment=u'#', index_col=0)
            HFB_Conv = (log_df['HFB_Conv'] == 'Yes')
            FAM_Conv = (log_df['FAM_Conv'] == 'Yes') if 'FAM_Conv' in log_df.columns else True
            Both_Conv = (HFB_Conv & FAM_Conv)
            if len(Both_Conv) > len(data):
                Both_Conv = Both_Conv[:len(data)]
            conv_info.extend(Both_Conv)
            nonConv_labels = [str(ii).zfill(6) for ii, conv in enumerate(Both_Conv) if not conv]
            if nonConv_labels:
                printf("Warning: PyNFAM run of", category, "has non-convergent point(s) with following label(s).")
                printf(*nonConv_labels)
            if 'Max_FAM_Si' in log_df:
                maxsi_info.extend(log_df['Max_FAM_Si'])
            else:
                maxsi_info.extend((np.nan,)*len(data))
            
        # save results into all_model_results and write it onto disk
        category_df = deepcopy(data)
        category_df['model_raw'] = model_results
        category_df['num_error_raw'] = num_errors
        df_func_residual(category_df, category) # apply func and calculate residuals
        category_df['convergence'] = conv_info
        category_df['max_fam_si'] = maxsi_info
        category_df['detail'] = output_dfs
        all_model_results.append(category_df)
        csv_cols = [col for col in category_df.columns if col != 'detail']
        category_df.to_csv(os.path.join(dir_path, result_filename + '_' + category + '.csv'), \
                           columns=csv_cols, index=False)
        category_df.to_pickle(os.path.join(dir_path, result_filename + '_' + category + '.pkl'))

        # clean replaced dir and files
        if clean_files:
            replaced_dir_path = os.path.join(scratch_dir_path, category+'_replaced')
            if os.path.isdir(replaced_dir_path):
                try:
                    move_or_del_file_dir(replaced_dir_path)
                except:
                    printf('Warning: Replaced dir', replaced_dir_path, 'is not fully cleaned!')
            try:
                files_rename_or_remove(dir_path, '.replaced', None)
            except:
                printf('Warning: Replaced files in', dir_path, 'are not fully cleaned!')
    # loop over categories ends here

    # finish check
    if check:
        if mu.do_mpi:
            comm.bcast(all_rerun_equal_zero, root=rank)
        return all_rerun_equal_zero

    # reorder and save results into output_data and also on disk
    output_data = pd.concat(all_model_results, sort=False)
    output_data.sort_index(inplace=True, ascending=True) # get back to the original order
    df_clean(output_data, 'exp_raw' in input_data) # drop unnecessary columns
    csv_cols = [col for col in output_data.columns if col != 'detail']
    output_data.to_csv(os.path.join(scratch_dir_path, result_filename +'.csv'), columns=csv_cols, index=False)
    output_data.to_pickle(os.path.join(scratch_dir_path, result_filename +'.pkl'))
    with open(os.path.join(scratch_dir_path, result_filename +'_other.pkl'), 'wb') as fileobj:
        pickle.dump(categories, fileobj)
        pickle.dump(all_rerun_equal_zero, fileobj)
        pickle.dump(start_from_backup, fileobj)
    if clean_files:
        try:
            # clean "result.pkl.replaced"
            for extension in ('.pkl', '.csv'):
                replaced_result_path = os.path.join(scratch_dir_path, result_filename + extension + '.replaced')
                if os.path.isfile(replaced_result_path):
                    os.remove(replaced_result_path)
        except:
            printf('Warning: Replaced result files are not fully cleaned!')
    
    # backup files
    if backup_option:
        overwrite_flag = False
        if (backup_option < 0) and (start_from_backup >= 0) and all_rerun_equal_zero:
            # backup to where we start from
            new_dir_name = str(start_from_backup)
            backup_num_limit(backup_max, backup_dir_path)
            overwrite_flag = True
        else:
            new_dir_name = backup_num_limit(backup_max, backup_dir_path)
        if abs(backup_option) == 1: # pkl file for check=True and replaced files are always ignored
            ignore_tuple = ('*.tar','*_check.pkl','*.replaced') # ignore detailed outputs in tars
        else:
            ignore_tuple = ('*_og.tar','*_check.pkl','*.replaced') # ignore old tars
        scratch_to_backup(os.path.join(backup_dir_path, new_dir_name), scratch_dir_path, ignore_tuple, \
                          categories, result_filename, overwrite_flag)

    # return
    if not return_df:
        output_data = output_data.to_dict('list')
    if mu.do_mpi:
        output_data = comm.bcast(output_data, root=rank)
    return output_data


def pynfam_fit_wrapper(
        input_data, categories=None, override_setts_fit={}, override_setts_cat={}, return_df=True, clean_files=False,
        exes_dir='./exes', scratch_dir='./scratch', backup_dir='./backup', backup_option=-1, backup_max=-1, start_from_backup=-1,
        nr_parallel_calcs=None, nthreads_per_calc=None, rerun_mode=None, ignore_hfb_nonconv=2, hfb_only=False, use_ratinterp=1,
        assert_consistency=False, comm=None, check=False, **kwargs):
    # This is the function open to outside for fitting
                   
    # Setup MPI variables
    if mu.do_mpi:
        if comm is None: # no MPI communicator input from outside wrapper
            comm = mu.MPI.COMM_WORLD
        assert(isinstance(comm, mu.MPI.Comm))
        rank, comm_size = mu.pynfam_mpi_traits(comm)
    else:
        rank, comm_size = 0, 1
    if mu.do_mpi and (comm_size == 1):
        mu.do_mpi = False
        print(u"*** Only one MPI task requested. Performing serial calculation. ***")

    if rank == 0:
        printf('Call of pynfam_fit_wrapper starts at', datetime.datetime.now())
    
    # parameters that cannot be changed by input parameters (for simplicity)
    # ignore_nonconv = 2
    # hfb_restart_file = None # HFB restart_file option
    # use_fam_storage = None # FAM use_fam_storage option
    # Ex_input = True # whether input energy is the excitation energy w.r.t. the daughter nucleus
    result_filename = 'result'

    # check consistency of parameters between processes
    if mu.do_mpi:
        assert_consistency = comm.bcast(assert_consistency, root=0)
    if assert_consistency and mu.do_mpi and comm_size > 1:
        flag = True
        my_params = locals()
        if rank == 0:
            all_params, all_input_data = zip(*comm.gather((None, None), root=0))
            for proc in range(1, comm_size):
                params = all_params[proc]
                input_data_recv = all_input_data[proc]
                for key, value in params.items():
                    flag = flag and (value == my_params[key])
                    if not flag: break
                if not flag: break
                if isinstance(input_data, pd.DataFrame) and isinstance(input_data_recv, pd.DataFrame):
                    flag = flag and (input_data.equals(input_data_recv))
                else:
                    flag = flag and (input_data == input_data_recv)
                if not flag: break
        else:
            to_send = {k: v for k, v in my_params.items() if k not in \
                    {'flag', 'assert_consistency', 'input_data', 'comm', 'rank', 'comm_size'}}
            comm.gather((to_send, input_data), root=0)
        flag = comm.bcast(flag, root=0)
        assert(flag)
    
    # main work
    if rank == 0:
        result = pynfam_fit_wrapper_root(input_data, categories, override_setts_fit, override_setts_cat, return_df, clean_files, \
                                         exes_dir, scratch_dir, backup_dir, backup_option, backup_max, start_from_backup, \
                                         nr_parallel_calcs, nthreads_per_calc, rerun_mode, ignore_hfb_nonconv, hfb_only, use_ratinterp, \
                                         comm, rank, check, result_filename, \
                                         **kwargs)
        printf('One run of pynfam_fit_wrapper finishes at', datetime.datetime.now())
    else:
        result = pynfam_fit_wrapper_nonroot(comm, 0)
    return result


class pynfam_residual(object):
    # class of residual functions for POUNDerS or other fitting / least_squares routines
    # Landau parameters and scaled isoscalar pairing should be the variables to fit

    result = None # result storage

    def __init__(
        self, input_data, categories=None, override_setts_fit={}, override_setts_cat={}, return_df=True, clean_files=False,
        exes_dir='./exes', scratch_dir='./scratch', backup_dir='./backup', backup_option=-1, backup_max=-1, start_from_backup=-1,
        nr_parallel_calcs=None, nthreads_per_calc=None, rerun_mode=None, ignore_hfb_nonconv=2, use_ratinterp=1, assert_consistency=False, comm=None, check=False,
        special_params_names=(), special_params_x0=(), special_params_reeval=lambda *args: args[0], special_params_bounds=(-np.inf, np.inf),
        special_params_xtol=1e-14, special_params_outfile=None, nonConv_penalty=None, posinf_replace=None, restart_in_seq=False, **kwargs):
        # initialize parameters that will not change during the fitting process
        self.args = deepcopy(locals())
        self.args.pop('self')
        self.args.update(self.args.pop('kwargs'))
        self.args['hfb_only'] = False
        self.__comm = self.args.pop('comm')

        # penalty for the non-converged point; if None, no penalty will be applied
        self.nonConv_penalty = self.args.pop('nonConv_penalty')

        # value to replace +inf
        self.posinf_replace = self.args.pop('posinf_replace')

        # special parameters whose change do not incur heavy calculations, default "GA"
        self.set_special_params(self.args.pop('special_params_names'), \
                                self.args.pop('special_params_x0'), \
                                self.args.pop('special_params_reeval'), \
                                self.args.pop('special_params_bounds'), \
                                self.args.pop('special_params_xtol'), \
                                self.args.pop('special_params_outfile'))
        
        # Setup MPI variables
        if mu.do_mpi:
            if self.comm is None: # no MPI communicator input from outside
                self.__comm = mu.MPI.COMM_WORLD
            assert(isinstance(self.comm, mu.MPI.Comm))
        self.__rank, _ = mu.pynfam_mpi_traits(self.comm)

        # prepare for restart of fit routine
        self.restart_in_seq = self.args.pop('restart_in_seq')
        if self.args['start_from_backup'] >= 0:
            # count starts from the given start_from_backup
            self.__count = self.args['start_from_backup']
        else:
            self.__count = 0
        self.__funflag = False
    
    @property
    def rank(self):
        # MPI process rank
        return self.__rank
    
    @property
    def comm(self):
        # MPI communicator
        return self.__comm
    
    @property
    def count(self):
        # count the number of calling "self.fun" to keep track of backups
        return self.__count
    
    @property
    def funflag(self):
        # flag to avoid direct calling of "call_pynfam" and "call_pynfam_extract_res"
        # when self.restart_in_seq is True.
        return self.__funflag

    def set_special_params(self, names, x0, routine, bounds=(-np.inf, np.inf), xtol=1e-14, outfile=None):
        # set names, starting point, bounds, evaluations function, tolerance and output file for special parameters
        # names and x0 must be an 1D arrays with the same length, and bounds must be a 2-tuple of 1D arrays. 
        assert((len(names) == len(x0)) and (len(bounds)==2))
        self.special_params_names = names
        self.special_params_vals = x0
        self.special_params_bounds = bounds
        self.special_params_reeval = routine
        self.special_params_xtol = xtol
        self.special_params_outfile = outfile
        if outfile: # print header line to outfile
            header_line = 'x    special_params    chi^2    ||g||_special    status'
            try:
                # if there exists old outfile with correct header, don't overwrite it
                with open(outfile, 'r') as fileobj:
                    first_line = fileobj.readline().rstrip()
                assert(first_line == header_line)
            except:
                # create or truncate outfile and write header line
                with open(outfile, 'w') as fileobj:
                    printf(header_line, file=fileobj)

    def call_pynfam(self, **fit_vars):
        # run pynfam_fit_wrapper with fitting variables and other parameters initialized in __init__
        if self.restart_in_seq and (not self.funflag):
            raise RuntimeError( \
                'pynfam_residual.call_pynfam must be called via pynfam_residual.fun with restart feature enabled!')
        allargs = deepcopy(self.args)
        allargs.update(fit_vars)
        self.result = pynfam_fit_wrapper(**allargs, comm=self.comm)
        return self.result

    @staticmethod
    def extract_cols(result, *cols):
        # extract and return columns of pynfam_fit_wrapper result in the form of a list of ndarrays
        r = [np.array(result[col]) for col in cols]
        if len(r) == 1:
            return r[0]
        else:
            return r

    def call_pynfam_extract_cols(self, *cols, **fit_vars):
        # run pynfam and extract specified columns / keys from the result
        if self.restart_in_seq:
            raise RuntimeError( \
                'pynfam_residual.call_pynfam_extract_cols cannot be called with restart feature enabled!')
        return self.extract_cols(self.call_pynfam(**fit_vars), *cols)

    @staticmethod
    def extract_residual(result, nonConv_penalty=None, posinf_replace=None):
        # extract and return model residuals in the form of ndarray, 
        # with penalty applied for non-convergent points and +inf replaced by a given value
        if 'norm_res' in result:
            r = np.array(result['norm_res'])
        else:
            printf('Warning: No normalized residual info obtained from pynfam_fit_wrapper! Fall back to residual.')
            r = np.array(result['residual'])
        if nonConv_penalty is not None:
            r = np.where(result['convergence'], r, (1.0+np.abs(nonConv_penalty))*r)
        else:
            r = np.where(result['convergence'], r, np.nan)
        if posinf_replace is not None:
            r = np.where(np.isposinf(r), posinf_replace, r)
        return r
    
    def call_pynfam_extract_res(self, **fit_vars):
        # run pynfam and extract residuals with penalty
        if self.restart_in_seq and (not self.funflag):
            raise RuntimeError( \
                'pynfam_residual.call_pynfam_extract_res must be called via pynfam_residual.fun with restart feature enabled!')
        return self.extract_residual(self.call_pynfam(**fit_vars), self.nonConv_penalty, self.posinf_replace)
    
    def _fun(self, xin, var_names):
        # auxiliary function that directly calls "call_pynfam_extract_res"
        assert(len(xin)==len(var_names))
        dict_fitvars = dict(zip(var_names, xin))
        return self.call_pynfam_extract_res(**dict_fitvars)

    def fun(self, xin, var_names):
        # function to be called by fitting / least_squares routines, with residuals as return
        # values and names of fitting variables are given by array "xin" and array "var_names", respectively
        if self.rank == 0:
            printf('Call pynfam_residual.fun with input args', xin, var_names)
        flag_modify = False # store whether start_from_backup is modified
        if self.restart_in_seq: # restart from existing solutions in backup
            assert(self.args.get('backup_option', -1) <= 0)
            backup_dir = self.args.get('backup_dir', './backup')
            if os.path.isdir(os.path.join(os.path.realpath(backup_dir), str(self.count))):
                # modify start_from_backup only when the backup of self.count exists
                flag_modify = True
                init_start_from_backup = self.args.get('start_from_backup', -1)
                self.args['start_from_backup'] = self.count
        self.__count += 1 # increase count
        self.__funflag = True # enable calling "call_pynfam_extract_res"
        residuals = self._fun(xin, var_names) # run pynfam and extract residuals
        self.__funflag = False # disable calling "call_pynfam_extract_res"
        if flag_modify: # change start_from_backup back
            self.args['start_from_backup'] = init_start_from_backup
        return residuals

    # some common fit choices
    def fun_g0p(self, xin):
        return self.fun(xin, ('g0p',))
    def fun_g0p_vpair0(self, xin):
        return self.fun(xin, ('g0p','vpair_t0_scaled'))
    def fun_g0p_g1p(self, xin):
        return self.fun(xin, ('g0p','g1p'))
    def fun_g0p_g1p_vpair0(self, xin):
        return self.fun(xin, ('g0p','g1p','vpair_t0_scaled'))
        
    def fun_special(self, xin, var_names):
        # function to be called by fitting / least_squares routines, with residuals as return
        # values and names of fitting variables are given by array "xin" and array "var_names", respectively
        # special parameters specified in self.special_params_* are not given via "xin" and "var_names"
        # instead, chi^2 is minimized separately in the direction of these special params, i.e. $\min_{special params} \chi^2$. 
        if self.rank == 0:
            printf('Call pynfam_residual.fun_special with input args', xin, var_names)
        if self.special_params_names:
            # avoid backup in the first run of pynfam
            backup_option = self.args.get('backup_option', -1)
            self.args['backup_option'] = 0
        # run pynfam
        residuals = self.fun(np.append(xin, self.special_params_vals), (*var_names, *self.special_params_names))
        if self.special_params_names: # deal with special params
            if self.rank == 0: # minimize chi^2 in the direction of special params
                def t(spe_params_vals): # auxiliary function for scipy.least_squares
                    self.result = self.special_params_reeval(self.result, **dict(zip(self.special_params_names, spe_params_vals)))
                    return self.extract_residual(self.result, self.nonConv_penalty, self.posinf_replace)
                printf('Scipy.least_squares in pynfam_residual.fun_special starts at', datetime.datetime.now())
                ls_res = least_squares(t, self.special_params_vals, bounds=self.special_params_bounds, \
                                       ftol=None, gtol=None, xtol=self.special_params_xtol)
                printf('Scipy.least_squares in pynfam_residual.fun_special finishes at', datetime.datetime.now())
                printf('Scipy least_squares summary for', self.special_params_names, ':')
                printf(ls_res)
                self.special_params_vals = ls_res.x
                if not ls_res.success:
                    printf('Warning: Scipy least_squares fails!')
                if self.special_params_outfile: # record to file
                    with open(self.special_params_outfile, 'a') as fileobj:
                        printf(np.array(xin).tolist(), ls_res.x.tolist(), 2*ls_res.cost, \
                               2*np.sqrt(np.sum(ls_res.grad**2)), ls_res.status, file=fileobj, sep='; ')
            if mu.do_mpi: # broadcast optimal special params
                self.special_params_vals = self.comm.bcast(self.special_params_vals, root=0)
            # backup arrangements
            init_start_from_backup = self.args.get('start_from_backup', -1) # store old start_from_backup
            self.args['start_from_backup'] = -1 # always start from scratch in the second pynfam run
            if self.rank == 0 and backup_option: # retrieve info for manual backup after the second run
                scratch_dir_path = os.path.realpath(self.args.get('scratch', './scratch'))
                result_filename = 'result'
                with open(os.path.join(scratch_dir_path, result_filename +'_other.pkl'), 'rb') as fileobj:
                    categories = pickle.load(fileobj)
                    all_rerun_equal_zero = pickle.load(fileobj)
                    start_from_backup = pickle.load(fileobj)
            # run pynfam with new values of special params so that we have corresponding files in the backup
            self.__funflag = True # enable calling "call_pynfam_extract_res"
            residuals = self._fun(np.append(xin, self.special_params_vals), (*var_names, *self.special_params_names))
            self.__funflag = False # disable calling "call_pynfam_extract_res"
            self.args['backup_option'] = backup_option # change backup_option back
            self.args['start_from_backup'] = init_start_from_backup # change start_from_backup back
            # backup manually
            if self.rank == 0 and backup_option:
                backup_max = self.args.get('backup_max', -1)
                backup_dir_path = os.path.realpath(self.args.get('backup_dir', './backup'))
                if not os.path.isdir(backup_dir_path):
                    os.makedirs(backup_dir_path)
                overwrite_flag = False
                if (backup_option < 0) and (start_from_backup >= 0) and all_rerun_equal_zero:
                    # backup to where we start from
                    new_dir_name = str(start_from_backup)
                    backup_num_limit(backup_max, backup_dir_path)
                    overwrite_flag = True
                else:
                    new_dir_name = backup_num_limit(backup_max, backup_dir_path)
                if abs(backup_option) == 1: # pkl file for check=True and replaced files are always ignored
                    ignore_tuple = ('*.tar','*_check.pkl','*.replaced') # ignore detailed outputs in tars
                else:
                    ignore_tuple = ('*_og.tar','*_check.pkl','*.replaced') # ignore old tars
                scratch_to_backup(os.path.join(backup_dir_path, new_dir_name), scratch_dir_path, ignore_tuple, \
                                  categories, result_filename, overwrite_flag)
        return residuals

    @staticmethod
    def static_find_minChiSqr_in_backups(var_names, backup_path='./backup', result_filename='result', find_in_range=None):
        # find the solution of minimum Chi square from existing solutions in the backup folder
        # var_names gives an array of strings with all the names of fitting parameters
        backup_path = os.path.realpath(backup_path)
        all_params = []
        all_chi_square = []
        all_res_vec = []
        all_nums = []
        for num in sorted([int(s) for s in os.listdir(backup_path) if s.isdigit()]):
            if (find_in_range is not None) and (num not in find_in_range):
                continue # only consider num in the given range
            inner_path = os.path.join(backup_path, str(num))
            if not os.path.isdir(inner_path): continue
            result_path = os.path.join(inner_path, result_filename + '.pkl')
            if not os.path.isfile(result_path): continue
            result_df = pd.read_pickle(result_path)
            prev_param_list = []
            for category in os.listdir(inner_path):
                category_path = os.path.join(inner_path, category)
                if not os.path.isdir(category_path): continue
                inputpkl_path = os.path.join(category_path, 'pynfam_inputs_'+category+'.pkl')
                if not os.path.isfile(inputpkl_path): continue
                with open(inputpkl_path, 'rb') as f:
                    pickle.load(f)
                    params_dict = pickle.load(f)
                param_list = []
                for param in var_names:
                    for inner_dict in params_dict.values():
                        if param in inner_dict.keys():
                            param_list.append(inner_dict[param])
                            break
                if prev_param_list: assert(param_list == prev_param_list)
                prev_param_list = param_list
            res_vec = pynfam_residual.extract_residual(result_df)
            if not np.isnan(res_vec).any():
                # avoid results with NaNs caused by runtime errors
                all_nums.append(num)
                all_params.append(np.array(param_list))
                all_res_vec.append(res_vec)
                all_chi_square.append(np.sum(all_res_vec[-1]**2))
        if all_params:
            zipped = zip(all_params, all_res_vec, all_chi_square, all_nums)
            min_zipped = min(zipped, key=lambda x: x[2])
            return min_zipped
        else:
            return (None,)*4

    def find_minChiSqr_in_backups(self, var_names, result_filename='result', find_in_range=None):
        return self.static_find_minChiSqr_in_backups(var_names, self.args.get('backup_dir', './backup'), result_filename, find_in_range)

    @staticmethod
    def vary_GA_only(df_in, GA, inplace=False):
        # update the DataFrame obtained from pynfam_fit_wrapper when only GA is varied.
        if isinstance(df_in, dict):
            # convert dict input to DataFrame
            save_dict = df_in
            df_in = pd.DataFrame.from_dict(df_in)
        else:
            save_dict = None
            assert(isinstance(df_in, pd.DataFrame))
        if inplace:
            df = df_in
        else:
            df = deepcopy(df_in)
        if 'model_raw' not in df: df['model_raw'] = df['model']
        if 'num_error_raw' not in df:
            if 'num_error' in df:
                df['num_error_raw'] = df['num_error']
            else:
                df['num_error'] = np.nan
                df['num_error_raw'] = np.nan
        exp_raw_in_input = ('exp_raw' in df)
        column_name = 'Half-Life(s)'
        for ii, (cat, detail_df) in enumerate(zip(df['category'], df['detail'])):
            if cat != 'HL': continue # only half life is impacted by GA
            if isinstance(detail_df, pd.DataFrame):
                temp_detail_df = detail_df if inplace else deepcopy(detail_df)
                detaildf_process_HL(temp_detail_df, GA=GA) # recalculate rates
                df['detail'].iat[ii] = temp_detail_df
                df['model_raw'].iat[ii] = temp_detail_df.at['Total', column_name]
                df['num_error_raw'].iat[ii] = temp_detail_df.at['Total', 'NumErr_'+column_name]
            else: # when there is no detail, fill with NaN
                df['model_raw'].iat[ii] = df['num_error_raw'].iat[ii] = np.nan
        df_func_residual(df, False) # apply func and calculate residuals
        df_clean(df, exp_raw_in_input) # clean NaN
        if save_dict is not None: # return dict if input is dict
            df = df.to_dict('list')
            if inplace:
                save_dict.clear()
                save_dict.update(df)
        if inplace:
            return None
        else:
            return df