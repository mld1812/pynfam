# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
from builtins   import zip
from builtins   import object
from builtins   import str
# -------------- Utilities -----------------
from copy        import deepcopy
from scipy.io    import FortranFile
from collections import defaultdict
import numpy  as np
import pandas as pd
import os
# ----------- Relative Imports -------------
from ..outputs.pnfam_parser import strengthOutParser
from ..fortran.pnfam_run import pnfamRun
from ..fortran.hfbtho_run import hfbthoRun
from ..utilities.hfb_utils import convString
from .contour            import famContour
from ..config import TMIN
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-07-26'

#===============================================================================#
#                             CLASS famStrength                                 #
#===============================================================================#
class famStrength(object):
    """
    An instance of the class famStrength contains all the info needed to write a .ctr
    binary file, namely an (operator + kval), a contour, and the nucleus (N,Z,A).
    It contains the methods to generate a list of fam objects to run, process their
    results, and write those results to a .out textfile and a .ctr binary file.

    These objects are manipulated by the shapeFactor class to compute shape
    factors and rates.

    Args:
        operator (str): Operator name as in the pnfam namelist.
        k (int): K-value of the operator as in the pnfam namelist.
        contour (famContour, str): String is assumed to be contour type, which
            instantiates a default contour of that type.
        nucleus (tuple): (N,Z,A) (default None). None is for if it's to be read
            from .ctr file.

    Attributes:
        op (str): Operator name as in the pnfam namelist.
        k (int): K-value of the operator as in the pnfam namelist.
        nucleus (tuple): (N,Z,A)
        contour (famContour): The complex energy contour along which strength is computed.
        str_df (dataframe): The strength function results (Re and Im columns are separated).
        meta_df (dataframe): The times and convergence of each fam point.
        _meta (dict): Internal initialization for the meta attribute.
    """

    def __init__(self, operator, k, contour, hfb=None):
        # Operator characteristics
        self.op    = operator
        self.k     = k

        # Nucleus and contour
        self.nucleus = None
        self.temperature = 0.0
        self.ft_active = False
        if hfb is not None:
            self.nucleus = hfb.nucleus
            self.temperature = hfb.temperature
            self.ft_active = hfb.ft_active
        if isinstance(contour, str):
            contour = famContour(contour)
        self.contour = contour

        # Strength result
        self.str_df = None
        self.meta_df = None
        self._meta = {u'Version':u'Unknown',
                     u'Interaction':u'Unknown',
                     u'Time': None,
                     u'Conv': None}


    # File names based off of operator
    @property
    def file_txt(self):
        """ file_txt (str): File name of textfile output.
        """
        return u"{:}.out".format(self.opname)
    @property
    def file_bin(self):
        """ file_bin (str): File name of binary output.
        """
        return u"{:}.out.ctr".format(self.opname)
    @property
    def file_tar(self):
        """ file_tar (str): File name of tarfile.
        """
        return u"{:}.tar".format(self.opname)

    # Different info/strings for the operator
    @property
    def opname(self):
        """ opname (str): Bare_operator + beta_type + kval.
        """
        return u"{:}K{:1d}".format(self.op, self.k)
    @property
    def bareop(self):
        """ bareop (str): Bare operator name without beta type.
        """
        return self.op[:-1]
    @property
    def beta(self):
        """ beta (str): Beta decay type.
        """
        return self.op[-1]
    @property
    def genopname(self):
        """ genopname (str): Bare_operator + _ + kval.
        """
        return u"{:}_K{:1d}".format(self.bareop, self.k)

    @property
    def xterms(self):
        """ xterms (list of str): The cross terms for this operator.
        """
        if self.str_df is None: return None
        return [str(c[3:-1]) for c in self.str_df.columns[2:] if u'Re' in c]
    @property
    def nxterms(self):
        """ nxterms (int): The number of cross terms for this operator.
        """
        if self.str_df is None: return None
        return len([c for c in self.str_df.columns[2:] if u'Re' in c])
    @property
    def meta(self):
        """ meta (dict): Meta data for the output header. Includes conv and time
        for every fam point from the meta_df attribute.
        """
        if self.meta_df is not None:
            self._meta[u'Conv']=self.concatFamConv(list(self.meta_df[u'Conv'].values))[0]
            self._meta[u'Time']=self.concatFamTime(list(self.meta_df[u'Time'].values))[0]
        return self._meta

    @property
    def ctr_param_df(self):
        """ ctr_param_df (dataframe): Contour parameterization (Theta, EQRPA)
        """
        headers1 = [u'Theta',u'Re(EQRPA)',u'Im(EQRPA)']
        df = pd.DataFrame()
        if self.contour.closed:
            df[headers1[0]] = self.contour.theta
        df[headers1[1]] = np.real(self.contour.ctr_z)
        df[headers1[2]] = np.imag(self.contour.ctr_z)
        return df

    @property
    def ctr_integ_df(self):
        """ ctr_ineg_df (dataframe): Contour integration parameters (GL_weights, dz/dt).
        """
        headers3 = [u'GL_Weights',u'Re(dzdt)',u'Im(dzdt)']
        df = pd.DataFrame()
        if self.contour.closed:
            df[headers3[0]] = self.contour.glwts
            df[headers3[1]] = np.real(self.contour.ctr_dzdt)
            df[headers3[2]] = np.imag(self.contour.ctr_dzdt)
        return df

    @property
    def cstr_df(self):
        """ cstr_df (complex dataframe): Reformat str_df with one column per
        strength+crossterm, with complex entries.
        """
        if self.str_df is None: return None
        str_cmp = pd.DataFrame()
        for i, col in enumerate(self.str_df):
            if i%2==0:
                str_real = self.str_df[col].values
            else:
                str_cmp[col[3:-1]] = str_real + self.str_df[col].values*1j
        return str_cmp

    @property
    def cstr_df_zeroed(self):
        """ cstr_df_zeroed (complex dataframe): cstr_df but any strength below cutoff
        is set to zero (including cross terms).
        """
        if self.str_df is None or self.contour.closed:
            return None
        cstrz = self.cstr_df.copy()
        cutoff = 1e-7 #self.contour.half_width
        cstrz.iloc[cstrz.index[np.imag(cstrz[u'Strength'].values)<cutoff],:] = 0
        return cstrz

    @property
    def use_FT_prefactor(self):
        """ use_FT_prefactor (bool): Indicates if strength contains the finite temp. prefactor
        """
        if self.ft_active:
            if self.temperature > TMIN or not self.contour.closed:
                return True
        return False

    #--------------------------------------------------------------
    def getFamList(self, paths, nml_params=None, subdir=u""):
        """
        Make a list of fam objects for each contour point and corresponding to
        this strength operator.

        Args:
            paths (pynfamPaths): pynfamPaths object to instantiate pnfamRun objects.
            nml_params (dict): Adjust pnfamRun default settings (default None).

        Returns:
            list of pnfamRun
        """

        # Make a representative object with operator, paths, label prefix, and updated nml
        fam_rep       = pnfamRun(paths, self.op, self.k)
        fam_rep.label = os.path.join(paths.rp(paths.fam_m),subdir,self.opname)
        if nml_params is not None:
            fam_rep.setNmlParam(nml_params)

        self.meta[u'Interaction'] = fam_rep.interaction

        # Assign contour points to copies of representative object
        fam_list     = []
        nr_compute   = self.contour.nr_compute
        computed_pts = self.contour.ctr_z[:nr_compute]
        for i, pt in enumerate(computed_pts):
            fam_pt = deepcopy(fam_rep)
            fam_pt.label = os.path.join(fam_rep.label, str(i).zfill(6))
            fam_pt.tmp_inputs.append(os.path.join(paths.hfb, hfbthoRun.file_bin))
            fam_pt.tmp_inputs.append(os.path.join(paths.hfb, hfbthoRun.file_nml))
            # Two body currents: If tbc file exists, use it as tmp input
            if fam_pt.nml[u'EXT_FIELD']['two_body_current_mode']:
                tbc_input = os.path.join(paths.hfb, u'{:}.tbc'.format(self.opname))
                if os.path.exists(tbc_input):
                    fam_pt.tmp_inputs.append(tbc_input)


            fam_pt.setCtrPoint(pt)
            fam_pt.ctr_params = {
                u'dzdt'  : self.contour.ctr_dzdt[i],
                u'theta' : self.contour.theta[i],
                u'glwt'  : self.contour.glwts[i]
                }
            fam_list.append(fam_pt)

        return fam_list

    #--------------------------------------------------------------
    def checkFamListCtr(self, fam_list):
        """
        Check a list of pnfamRun objects correspond to the same contour as this
        strength object. Attempt to fix before raising error.

        Args:
            fam_list (list of pnfamRun): The list to be checked.

        Raises:
            RuntimeError
        """

        # NB!: If we read f.omega from existing namelist, the float might
        # not exactly match that in ctr_z, so check within 1e-10...
        err = None
        if any(f.opname != self.opname for f in fam_list):
            err = u"A fam object operator does not match strength object operator."
        if len(fam_list) != self.contour.nr_compute:
            err = u"len(fam_list) != nr_compute."
        if len(fam_list) != self.contour.nr_points:
            if len(fam_list) != ((self.contour.nr_points + 1)//2):
                err = u"Invalid len(fam_list), must be nr_points or (nr_points+1)/2."
        if any([abs(np.real(f.omega-self.contour.ctr_z[i]))>1e-10 or\
                abs(np.imag(f.omega-self.contour.ctr_z[i]))>1e-10 \
                for i, f in enumerate(fam_list)]):
            # Try to fix by sorting before raising error
            try:
                fam_map    = {"{:.10e}".format(f.omega) : f for f in fam_list}
                fam_sorted = [fam_map["{:.10e}".format(z)] for z in self.contour.ctr_z[:self.contour.nr_compute]]
                if any([abs(np.real(f.omega-self.contour.ctr_z[i]))>1e-10 or\
                        abs(np.imag(f.omega-self.contour.ctr_z[i]))>1e-10 \
                        for i, f in enumerate(fam_list)]):
                    raise RuntimeError
                fam_list = fam_sorted
            except:
                err = u"fam omega != contour ctr_z."
        if err is not None:
            raise RuntimeError(u"fam data does not match contour data.\n",err)

    #--------------------------------------------------------------
    def concatFamStr(self, fam_list):
        """
        Combine data from a list of pnfamRun objects into strength dataframe.
        If contour is closed and nr_compute = (nr_points+1)/2, fill in missing
        points with symmetry S(w) = S*(w*).

        Args:
            fam_list (list of pnfamRun): List to be compiled.

        Returns:
            DataFrame: Contains just strength and cross terms, no contour data.
        """

        self.checkFamListCtr(fam_list)

        # DF with just strength and xterms (partial if ctr is closed)
        op_str_dict = defaultdict(list)
        for fam_pt in fam_list:
            for str_pt in fam_pt.soln_dict[u'Strength']:
                op_str_dict[str_pt[u'label']].append(str_pt[u'val'])
        op_str_df = pd.DataFrame(op_str_dict)
        header_order = [s[u'label'] for s in fam_list[0].soln_dict[u'Strength']]

        # Fill in missing strength points with symmetry s(w) = s*(w*)
        if len(fam_list) != self.contour.nr_points:
            endpoints  = self.contour.ctr_z[0] == self.contour.ctr_z[-1]
            nr_compute = self.contour.nr_compute
            nr_points  = self.contour.nr_points

            # flip row order and copy
            sym_str_df = op_str_df.iloc[::-1].copy()

            # avoid double counting the last point if nr_points is odd
            if nr_compute*2 != nr_points:
                sym_str_df.drop(sym_str_df.head(1).index, inplace=True)

            # Conjugate (s(w) = s*(w*))
            for h in header_order:
                if u'Im' in h: sym_str_df[h] = -sym_str_df[h].values

            # Fudge endpoints if its the same point (e.g. evenly spaced grid)
            if endpoints and nr_points != 1:
                sym_str_df.iloc[-1] = op_str_df.iloc[0].values

            op_str_df = pd.concat([op_str_df, sym_str_df], axis=0, ignore_index=True)

        # Give columns same order as in original pnfam output
        op_str_df = op_str_df[header_order]

        # Include FT prefactor (DONT for closed contours with T <= TMIN)
        if self.ft_active:
            if self.temperature > TMIN or not self.contour.closed:
                temp = self.temperature
                ctr_z = self.contour.ctr_z
                # Pretend open contours are on the real axis
                if not self.contour.closed:
                    pf = 1.0/(1.0 - np.exp(-np.real(ctr_z)/temp))
                else:
                    pf = 1.0/(1.0 - np.exp(-ctr_z/temp))
                pf[np.logical_not(np.isfinite(pf))] = 0 # Is this desired?

                # Assume headers are Re,Im,Re,Im,... but don't Re,Im vs Im,Re
                for i in range(len(header_order)//2):
                    h1 = header_order[i*2]
                    h2 = header_order[i*2+1]
                    s1 = op_str_df[h1]
                    s2 = op_str_df[h2]
                    if u'Im' in h1:
                        s = (s2 + s1*1j)*pf
                        op_str_df[h1] = np.imag(s)
                        op_str_df[h2] = np.real(s)
                    else:
                        s = (s1 + s2*1j)*pf
                        op_str_df[h2] = np.imag(s)
                        op_str_df[h1] = np.real(s)

        return op_str_df

    #--------------------------------------------------------------
    def concatFamTime(self, fam_list):
        """
        Get the time data from a list of pnfamRun solutions.

        Args:
            fam_list (list of pnfamRun or float): List to be compiled.

        Returns:
            float: The total time.
            ndarray: The times for each pnfamRun.
        """
        if all(isinstance(f,pnfamRun) for f in fam_list):
            time_list = [f.soln_dict[u'Time'] for f in fam_list]
        else:
            time_list = fam_list
        times = np.array(time_list)
        total_time = np.sum(times) # pnfam2 times are in minutes
        return total_time, times

    #--------------------------------------------------------------
    def concatFamConv(self, fam_list):
        """
        Get the convergence data from a list of pnfamRun solutions.

        Args:
            fam_list (list of pnfamRun or str): List to be compiled.

        Returns:
            str: The net convergence.
            ndarray: The convergence for each pnfamRun.
        """
        if all(isinstance(f,pnfamRun) for f in fam_list):
            conv_list = [f.soln_dict[u'Conv'] for f in fam_list]
        else:
            conv_list = fam_list

        conv = convString(conv_list)

        return conv, conv_list

    #--------------------------------------------------------------
    def getMeta(self, data):
        """
        Get the meta data for output header from pnfam solutions. Store
        in meta_df and meta attributes.

        Args:
            data (list of pnfamRun or str): List to be compiled, or path
                to existing OP.out file created by pynfam.

        Returns:
            str: The net convergence.
            ndarray: The convergence for each pnfamRun.

        Raises:
            TypeError
            IOError
        """
        err = False
        if isinstance(data, list):
        # Assume list of fam objects
            conv_list = self.concatFamConv(data)[1]
            time_list = self.concatFamTime(data)[1]
            version = data[0].soln_dict[u'Version']
        elif isinstance(data, str):
        # Assume path to op.out file and parse
            sop = strengthOutParser(data)
            conv_list = sop.getAllConv()
            time_list = sop.getAllTime()
            version = sop.getVersion()
            err = sop.err
        else:
            raise TypeError(u"Invalid input type for method famStrength.getMeta")

        # Set the (potentially error) values
        self.meta[u'Version'] = version
        self.meta_df = pd.DataFrame({u'Time':time_list, u'Conv':conv_list})
        if err: raise IOError(sop.err_msg)

    #--------------------------------------------------------------
    def concatFamData(self, fam_list):
        """
        Compile fam strength and meta data from list of pnfamRun objects.

        Args:
            fam_list (list of pnfamRun): List to be compiled.
        """
        self.str_df = self.concatFamStr(fam_list)
        self.getMeta(fam_list)

    #--------------------------------------------------------------
    def writeStrengthOut(self, dest=u'./', fname=None):
        """
        Write the OP.out pynfam summary file.

        Args:
            dest (str): File destination (default './').
            fname (str): File name (default None --> self.opname+'.out')
        """

        if self.str_df is None:
            raise RuntimeError(u"Strength dataframe has not been populated. Cannot write to .out")

        df_list = [self.ctr_param_df, self.str_df, self.ctr_integ_df, self.meta_df]
        out_df = pd.concat(df_list, axis=1)

        if fname is None: fname = self.file_txt
        pd_string = out_df.to_string(header=True, index=True, col_space=3,
                            float_format=lambda x: u'{:25.16e}'.format(x))

        header = [u"# pnFAM code version:   {:}\n".format(self.meta[u'Version']),
                  u"# Total run time:       {:<.6f} mins\n".format(self.meta[u'Time']),
                  u"# All points converged: {:}\n".format(self.meta[u'Conv']),
                  u"# Residual interaction: {:}\n".format(self.meta[u'Interaction']),
                  u"# Operator:             {:} with K={:d}\n".format(self.op, self.k),
                  u"# Contour:              {:}\n".format(self.contour.name_and_int),
                  u"# Temperature:          {:<.6f} (Prefactor={:})\n".format(self.temperature, self.use_FT_prefactor),
                  u"#\n"]

        with open(os.path.join(dest, fname), u'w') as str_file:
            str_file.writelines(header)
            str_file.write(pd_string+u'\n')

    #--------------------------------------------------------------
    def writeCtrBinary(self, dest=u'./', fname=None):
        """
        Write an OP.ctr binary file like the one generated by pnfam contour
        mode that is compatible with the fortran betadecay.x.

        Args:
            dest (str): File destination (default './').
            fname (str): File name (default None --> self.opname+'.out.ctr')

        Raises:
            RuntimeError
        """

        # Types - hopefully identical to those used in fortran
        int_type = np.int32
        dbl_type = np.float_
        cmp_type = np.complex_
        cha_type = np.bytes_

        version = 3

        if self.str_df is None:
            raise RuntimeError(u"Strength dataframe has not been populated. Cannot write to .ctr")

        # Use the scipy FortranFile object to write and unformatted binary
        if fname is None: fname = self.file_bin
        ctr_bin =FortranFile(os.path.join(dest, fname), mode='w')

        # integer               file_version
        #-------------------------------------------------------------------------
        ctr_bin.write_record(np.array([version], dtype=int_type))
        # integer               Z(parent), A
        #-------------------------------------------------------------------------
        ctr_bin.write_record(np.array(\
                [self.nucleus[1], self.nucleus[2]], dtype=int_type))
        # character(80)         Operator
        #-------------------------------------------------------------------------
        formatted_opname = self.bareop.encode() + (80-len(self.bareop))*b' '
        ctr_bin.write_record(np.array([formatted_opname], dtype=cha_type))
        # integer               K, nr_points, nxterms
        #-------------------------------------------------------------------------
        ctr_bin.write_record(np.array([self.k, self.contour.nr_points, self.nxterms],
            dtype=int_type))
        # character(80) array   xterm_labels(:)
        #-------------------------------------------------------------------------
        if self.nxterms > 0:
            # xterms are already in bareop form so we're okay there
            formatted_xterms = [xt.encode()+(80-len(xt))*b' ' for xt in self.xterms]
        else:
            formatted_xterms = [80*b' ']
        xlabels = np.array(formatted_xterms, dtype=cha_type)
        ctr_bin.write_record(xlabels)
        # double array          ctr_t(nr_points)
        #-------------------------------------------------------------------------
        ctr_bin.write_record(self.contour.theta.astype(dbl_type))
        # dcomplex array        ctr_dzdt(nr_points)
        #-------------------------------------------------------------------------
        ctr_bin.write_record(self.contour.ctr_dzdt.astype(cmp_type))
        # dcomplex array        ctr_z(nr_points)
        #-------------------------------------------------------------------------
        ctr_bin.write_record(self.contour.ctr_z.astype(cmp_type))
        # dcomplex array        ctr_strength(nr_points, nxterms)
        #-------------------------------------------------------------------------
        # python stores row ordered but fortran is column ordered so write transpose
        ctr_bin.write_record((self.cstr_df.values.T).astype(cmp_type))
        # integer               use_gauleg_contour (0=no, 1=yes)
        #-------------------------------------------------------------------------
        ctr_bin.write_record(np.array([int(self.contour.use_gauleg)],dtype=int_type))
        # double array          ctr_gaugleg_weights(nr_points)
        #-------------------------------------------------------------------------
        ctr_bin.write_record(self.contour.glwts.astype(dbl_type))

        ctr_bin.close()

    #--------------------------------------------------------------
    def readCtrBinary(self, src=u'./', fname=None):
        """
        Read a .out.ctr binary file and use what is read to populate instance attributes.

        Args:
            src (str): File location (default './').
            fname (str): File name (default None --> self.opname+'.out')
        """

        int_type = np.int32
        dbl_type = np.float_
        cmp_type = np.complex_
        cha_type = np.bytes_

        if fname is None: fname  = self.file_bin
        f2read = os.path.join(src,fname)
        if not os.path.exists(f2read): raise IOError(u"Binary file not found.")

        # Read in the raw data (In pynfam V2.0 all strings are len=80, but in V1.0
        # if len(xterms)=0 we wrote str('') with len=1. When fortranfile raises error
        # it doesn't deal well with partially read contents, so just re-read the
        # whole file from scratch.
        for char in [80,1]:
            f = FortranFile(f2read, mode=u'r')
            version    = f.read_ints(dtype=int_type)[0]
            nucleus    = f.read_ints(dtype=int_type)
            op         = f.read_record(dtype=(cha_type,80))[0]
            k_npts_nxts= f.read_ints(dtype=int_type)
            try:
                xterms = f.read_record(dtype=(cha_type,char))
                break
            except ValueError as e:
                f.close()
        theta      = f.read_reals(dtype=dbl_type)
        ctr_dzdt   = f.read_record(dtype=cmp_type)
        ctr_z      = f.read_record(dtype=cmp_type)
        strength_raw = f.read_record(dtype=cmp_type)
        use_gauleg = f.read_ints(dtype=int_type)[0]
        glwts      = f.read_reals(dtype=dbl_type)
        f.close()

        # Format the data
        nuc       = (nucleus[1] - nucleus[0], nucleus[0], nucleus[1])
        xterms    = [xt.decode() for xt in xterms]
        xterms    = [xt.strip() for xt in xterms if not xt.isspace()]
        strength  = strength_raw.reshape(len(xterms)+1, len(theta)).T
        k         = k_npts_nxts[0]
        nr_points = k_npts_nxts[1]
        nxterms   = k_npts_nxts[2]
        bareop    = op.decode().strip()
        str_df  = pd.DataFrame(); headers = [u'Strength']+xterms
        for h, col in zip(headers, strength.T):
            str_df[u'Re({:})'.format(h)] = np.real(col)
            str_df[u'Im({:})'.format(h)] = np.imag(col)

        # Check the data is for the correct operator (no way to check if beta type is correct...)
        if bareop != self.bareop or k != self.k:
            raise IOError(u"File contents do not match requested operator.\n"+\
                    u"Expected: {:}, {:}, Read: {:}, {:}.".format(self.bareop, self.k, bareop, k))

        # Artificially define the contour properties from what is read
        # (we assume closed and quad properties based on contour_type provided)
        closed = self.contour.closed
        quad   = self.contour.quadrature
        if any(t != 0.0 for t in theta) and not self.contour.closed:
            print(u"WARNING: Specified contour is not closed, but data appears to be closed.")
            print(u"         Manually changing to closed contour.")
            closed = True; quad   = u'GAUSS'
        if closed:
            nr_compute = ((nr_points+1)//2)
        else:
            nr_compute = nr_points
        if any(z != np.imag(ctr_z[0]) for z in np.imag(ctr_z)):
            half_width = None
        else:
            half_width = np.imag(ctr_z[0])
        ctr_data = {
            u'nr_points'   : nr_points,
            u'nr_compute'  : nr_compute,
            u'use_gl_ctr'  : use_gauleg,
            u'ctr_z'       : ctr_z,
            u'ctr_dzdt'    : ctr_dzdt,
            u'theta'       : theta,
            u'glwts'       : glwts,
            u'half_width'  : half_width,
            u'quad'        : quad,
            u'closed'      : closed
            }

        # Populate attributes of the famStrength instance object
        self.str_df  = str_df # takes care of xterms, nxterms, and cstr_df via properties
        self.nucleus = nuc
        self.contour._ctr_data = ctr_data
