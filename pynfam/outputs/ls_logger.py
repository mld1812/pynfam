# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
from builtins   import object
# -------------- Utilities -----------------
import pandas as pd
import numpy as np
import os
import copy
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-07-26'

#===============================================================================#
#                             CLASS lsLogger                                    #
#===============================================================================#
class lsLogger(object):
    """
    An lsLogger contains a logfile name and the methods to format data, and read
    from and write to the logfile using pandas.

    Args:
      filename (str): The logfile name.
    """

    order_nuc = [u'Label', u'Z', u'Z-Blk', u'N', u'N-Blk']
    """ order_nuc (list of str): The order of columns in the logfile.
    """

    def __init__(self, filename=u'logfile.dat'):
        self.filename = filename

    #-----------------------------------------------------------------------
    def writeLog(self, fhandle, frame, output_dir=None, float_fmt=u'{:.6f}'):
        """
        Write a logfile for the provided dataframe.

        Args:
            fhandle (file): File handle object of the logfile.
            frame (dataframe): The data to be written to file.
            output_dir (str): Destination to which file is being written. Will be
                printed as a comment at the top of the logfile.
            float_fmt (str): Format string for floats in the data.
        """

        # Write the header line with the source directory
        if output_dir is not None:
            output_dir_str = u"# Logfile for contents contained in: " + output_dir
            fhandle.write(output_dir_str + u'\n')

        # Write the data titles/values
        if float_fmt is not None:
            float_format = lambda x: float_fmt.format(x)
        else:
            float_format = None

        # Replace None with np.nan
        frame.fillna(value=np.nan, inplace=True)
        # Pandas automatically truncates strings at 50 chars, as of now
        # the way to fix this is through the display options.
        pd.set_option('display.max_colwidth', 100)

        pd_str = frame.to_string(header=True, index=True, col_space=3,
                                 float_format=float_format)
        fhandle.write(pd_str+u'\n')

    #-----------------------------------------------------------------------
    def readLog(self, file_path=u'./', header=0, comment=u'#'):
        """
        Read in a logfile.

        Args:
            file_path (str): Path to logfile (default './')
            header (int): Line of the logfile to be used as header (Default 0).
            comment (str): Comment symbol in logfile (default '#').

        Returns:
            DataFrame
        """

        file2open = os.path.join(file_path, self.filename)
        df = pd.read_csv(file2open, delim_whitespace=True, header=header, comment=comment)
        return df

    #-----------------------------------------------------------------------
    def format2Frame(self, data):
        """
        Given dict or dataframe, return dataframe.

        Args:
            data (dict, dataframe): The data to be formatted.

        Returns:
            DataFrame
        """

        if isinstance(data, dict):
            dict4frame = copy.copy(data)
            for key,val in list(dict4frame.items()):
                if not isinstance(val, list):
                    dict4frame[key] = [val]
            data = pd.DataFrame(dict4frame)
        return data

    #-----------------------------------------------------------------------
    def orderedKeys(self):
        """
        Return list of keys in desired order.

        Returns:
            list
        """

        return lsLogger.order_nuc

    #-----------------------------------------------------------------------
    def orderedFrame(self, frame, order):
        """
        Reorder dataframe columns to begin with those in 'nuc_order', with
        extra column labels not in 'nuc_order' moved to the end.

        Args:
            frame (dataframe): The data to be re-ordered.
            order (list of str): The desired order of the columns.

        Returns:
            DataFrame
        """

        # Order the keys we have, allowing for some to be missing
        allkeys = list(frame.columns)
        ordered = [key for key in order if key in allkeys]
        # Get the remaining keys
        unordered = [key for key in allkeys if key not in ordered]

        ordered_keys = ordered + unordered

        # Remove debug entries
        if u'Debug' in allkeys:
            frame.drop([u'Debug'],inplace=True,axis=1)

        return frame[ordered_keys]

    #-----------------------------------------------------------------------
    def formatLog(self, data):
        """
        Convert data to dataframe and order the columns.

        Args:
            Data (dict, dataframe): The data to be formatted.

        Returns:
            DataFrame
        """
        frame = self.format2Frame(data)
        order = self.orderedKeys()
        return self.orderedFrame(frame, order)

    #-----------------------------------------------------------------------
    def quickWrite(self, logdata, dest=u'./', fmt=True, **kwargs):
        """
        Format then write a dictionary or dataframe to a logfile.

        Args:
            logdata (dict, dataframe): The data to be written.
            dest (str): Destination path.
            fmt (bool): If False, don't format the data.
            **kwargs: Keyword args for writeLog method.
        """
        file2open = os.path.join(dest, self.filename)
        if fmt:
            logdata = self.formatLog(logdata)
        with open(file2open, u'w') as f:
            self.writeLog(f, logdata, **kwargs)

#===============================================================================#
#                             CLASS hfbLogger                                   #
#===============================================================================#
class hfbLogger(lsLogger):
    """
    An hfbLogger is a subclass of the lsLogger class with extended methods for
    formatting the HFB logfile information.
    """

    order_hfb = [u'HFB_Conv', u'HFB_Time', u'Time', u'Energy', u'Def', u'KP']
    """ order_hfb (list of str): The order of HFB data in the logfile.
    """
    order_hfbbeta = [u'HFB_Qval',u'EQRPA_max',u'E_gs']
    """ order_hfbbeta (list of str): The order of HFB data requiring beta type.
    """

    #-----------------------------------------------------------------------
    def orderedKeys(self, frame):
        """
        Extend the lsLogger method to include the order of the HFB data.

        Args:
            frame (dataframe): The data to be written.

        Returns:
            list of str: The order of the columns.
        """

        order =  list(super(hfbLogger, self).orderedKeys())
        order += list(hfbLogger.order_hfb)
        # QP data only populated if hfbthoRun.beta is not None
        if u'HFB_Qval' in list(frame.columns):
            order += list(hfbLogger.order_hfbbeta)
        return order

    #-----------------------------------------------------------------------
    def formatLog(self, data):
        """
        Override the lsLogger method to format the HFB data.

        Args:
            data (dict, dataframe): The data to be formatted.

        Returns:
            DataFrame
        """

        frame = self.format2Frame(data)
        # Favor Total HFB Time in the log over time for a single run
        keys = list(frame.columns)
        if u'Total_Time' in keys and u'Time' in keys:
            frame.drop(columns=[u'Time'], inplace=True)
        frame.rename(columns={u'Conv':u'HFB_Conv', u'Total_Time':u'HFB_Time'}, inplace=True)
        order = self.orderedKeys(frame)
        return self.orderedFrame(frame, order)

#===============================================================================#
#                             CLASS betaLogger                                  #
#===============================================================================#
class betaLogger(hfbLogger):
    """
    A betaLogger is a subclass of the hfbLogger class with extended methods for
    formatting the beta decay logfile information.
    """

    order_beta = [u'FAM_Conv', u'FAM_Time', u'FAM_Qval', u'Total-HL', u'%FF']
    """ order_beta (list of str): The order of beta decay data in the logfile.
    """

    def orderedKeys(self, frame):
        """
        Extend the lsLogger method to include the order of the HFB data.

        Args:
            frame (dataframe): The data to be written.

        Returns:
            list of str: The order of the columns.
        """

        order = super(betaLogger, self).orderedKeys(frame)
        order += betaLogger.order_beta
        return order