# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
from builtins   import zip
# -------------- Utilities -----------------
import pandas as pd
import numpy as np
# ----------- Relative Imports -------------
from .parser import parser
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2020-08-18'

#===============================================================================#
#                             CLASS pnfamParser                                 #
#===============================================================================#
class pnfamParser(parser):
    """
    A pnfamParser is a subclass of the parser class that contains a
    pnfam output file and the methods to parse it for values.

    Args:
        pnfam_output (str, list): A str is assumed to be path to output file, which
            is then read in and converted to a list of strings. Newline symbols
            are stripped.

    Attributes:
        keys (dict): A collection of strings to search for in the output.
    """
    def __init__(self, pnfam_output):
        parser.__init__(self, pnfam_output)
        self.keys = {\
                     u'time'       : u'Total CPU time',
                     u'conv_yes'   : u'iteration converged',
                     u'conv_no'    : u'iteration interrupted',
                     u"version"    : u"Version:",
                     u"interaction": u"Residual interaction",
                     u"operator"   : u"Operator:",
                     u"half_width" : u"Gamma",
                     u"strength"   : u"Re(EQRPA)"}

    #-----------------------------------------------------------------------
    def getTime(self):
        """
        Parse for the time.

        Returns:
            float
        """
        try:
            time_line = self.getLineIndices(self.keys[u'time'])[0]
            time = self.getNumbers(self.output[time_line], -2, err_val=self.float_err)
        except Exception:
            time = self.float_err
            self.flagError(u"Error parsing time from FAM output.")
        return time

    #-----------------------------------------------------------------------
    def getConv(self):
        """
        Parse for the convergence.

        Returns:
            str
        """
        conv_yes = self.getLineIndices(self.keys[u'conv_yes'])
        conv_no  = self.getLineIndices(self.keys[u'conv_no'])
        if conv_yes:
            conv = u'Yes'
        elif conv_no:
            conv = u'No'
        else:
            conv = self.str_err
            self.flagError(u"No convergence lines found in FAM output.")
        return conv

    #-----------------------------------------------------------------------
    def getVersion(self):
        """
        Parse for the pnfam version.

        Returns:
            str
        """
        try:
            l_v = self.getLineIndices(self.keys[u'version'])[0]
            version = self.output[l_v].split()[-1].strip()
        except Exception:
            version = self.str_err
            self.flagError(u"Error parsing version from FAM output.")
        return version

    #-----------------------------------------------------------------------
    def getStrength(self):
        """
        Parse for the strength and cross term results.

        Returns:
            list of dict: Of form {'label: str, 'val': float} for each column.
        """
        try:
            df = pd.read_csv(self.src, delim_whitespace=True, header=0, comment=u'#', index_col=0)
            df.drop(df.index[0], inplace=True) # Drop the energy
            # Labels = Re(Strength) Im(Strength) Re(Xterm) Im(Xterm)...
            labels   = [ir+lab+")" for lab in list(df.index) for ir in [u"Re(", u"Im("]]
            str_vals = [None]*(2*len(df.index))
            str_vals[::2] = df[u"Real"].values
            str_vals[1::2] = df[u"Imag"].values
            if any(v is None for v in str_vals): raise Exception
            strength = [{u'label':h,u'val':float(d)} for h, d in zip(labels,str_vals)]
        except Exception:
            strength = 2*[{u'label':self.str_err, u'val':self.float_err}]
            self.flagError(u"Error parsing strength from FAM output.")
        return strength

    #-----------------------------------------------------------------------
    def getXterms(self):
        """
        Parse for the strength and cross term results.

        Returns:
            list of dict: Of form {'label: str, 'val': float} for each column.
        """
        try:
            df = pd.read_csv(self.src, delim_whitespace=True, header=0, comment=u'#', index_col=0)
            df.drop(df.index[0], inplace=True) # Drop the energy
            xterms   = list(df.index)[1:]
        except Exception:
            xterms = []
            self.flagError(u"Error parsing cross terms from FAM output.")
        return xterms

#===============================================================================#
#                             CLASS strengthOutParser                           #
#===============================================================================#
class strengthOutParser(parser):
    """
    A strengthOutParser is a subclass of the parser class that contains a pynfam
    compiled strength output file (OP.out) and the methods to parse it for values.

    Args:
        pnfam_output (str, list): A str is assumed to be path to output file, which
            is then read in and converted to a list of strings. Newline symbols
            are stripped.

    Attributes:
        keys (dict): A collection of strings to search for in the output.
    """
    def __init__(self, pnfam_output):
        parser.__init__(self, pnfam_output)
        self.keys = {u"version": u"pnFAM code version",
                     u"time"   : u"Total run time",
                     u"conv"   : u"All points converged"}

    #-----------------------------------------------------------------------
    def getVersion(self):
        """
        Parse for the pnfam version.

        Returns:
            str
        """
        try:
            l_v = self.getLineIndices(self.keys[u'version'])[0]
            version = self.output[l_v].split(u'version')[-1].strip()
        except Exception:
            version = self.str_err
            self.flagError(u"Error parsing version from strength output.")
        return version

    #-----------------------------------------------------------------------
    def getTime(self):
        """
        Parse for the total fam time.

        Returns:
            float
        """
        try:
            l_t = self.getLineIndices(self.keys[u'time'])[0]
            time = self.getNumbers(self.output[l_t], -2, err_val=self.float_err)
        except Exception:
            time = self.float_err
            self.flagError(u"Error parsing total time from strength output.")
        return time

    #-----------------------------------------------------------------------
    def getConv(self):
        """
        Parse for the total fam convergence.

        Returns:
            str
        """
        try:
            l_c = self.getLineIndices(self.keys[u'conv'])[0]
            conv = self.output[l_c].split()[-1].strip()
        except Exception:
            conv = self.str_err
            self.flagError(u"Error parsing total convergence from strength output.")
        return conv

    #-----------------------------------------------------------------------
    def getAllConv(self):
        """
        Parse for the entire convergence column.

        Returns:
            ndarray
        """
        try:
            df = pd.read_csv(self.src, delim_whitespace=True, header=0, comment=u'#')
            conv_list = df[u'Conv'].dropna().values
        except Exception:
            conv_list = np.array([self.str_err])
            self.flagError(u"Error parsing convergence list from strength output.")
        return conv_list

    #-----------------------------------------------------------------------
    def getAllTime(self):
        """
        Parse for the entire time column.

        Returns:
            ndarray
        """
        try:
            df = pd.read_csv(self.src, delim_whitespace=True, header=0, comment=u'#')
            time_list = df[u'Time'].dropna().values
        except Exception:
            time_list = np.array([self.float_err])
            self.flagError(u"Error parsing time list from strength output.")
        return time_list
