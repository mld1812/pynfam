# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
from builtins   import filter
from builtins   import object
from builtins   import str
# -------------- Utilities -----------------
import numpy as np
import re
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-07-26'

#===============================================================================#
#                             CLASS parser                                      #
#===============================================================================#
class parser(object):
    """
    An instance of the class parser contains the contents of an output file
    (as a list of strings) and the methods to parse it for values.

    Args:
        output_src (str, list): A str is assumed to be path to output file, which
            is then read in and converted to a list of strings. Newline symbols
            are stripped.

    Attributes:
        float_err (float, nan): Default error value for floats.
        int_err (int): Default error value for ints.
        str_err (str): Default error value for strings.
        err (bool): Were any errors encountered by the parser?
        err_msg (str): Error string if errors were encountered.
    """

    # [-+]?  = optional sign
    # \ *    = 0 or more spaces
    # [0-9]* = 0 or more digits
    # \.?    = optional dot ('\.' searches for a dot while '.' is a wildcard)
    # [0-9]+ = 1 or more digits
    # (?:)   = non-capturing group
    # ([eE]...)? = optional (e or E + optional sign + 1 or more digits)
    reg_ex =r'[-+]?\ *[0-9]*\.?[0-9]+(?:[eE]\ *[-+]?\ *[0-9]+)?'
    """ The regular expression used to parse for numbers.
    """

    def __init__(self, output_src):
        self.output    = []
        self.src       = None
        self.updateOutput(output_src)
        self.float_err = np.nan
        self.int_err   = 999
        self.str_err   = u'Error'
        self.err       = False
        self.err_msg   = u''

    #-----------------------------------------------------------------------
    def updateOutput(self, output_src):
        """
        Set the parser output attribute.

        * If output=str, assume it's the path to the output file and read.
        * If output=list, assume it's a list of strings (1 per line).
        * Invalid output reverts to empty list, which should return error values
          when any method is invoked.

        Args:
            output_src (str, list): The output.

        Raises:
            ValueError
        """
        if output_src is None:
            return
        elif isinstance(output_src, str):
            self.src = output_src
            try:
                with open(output_src, u'r') as out_file:
                    out_data = out_file.readlines()
                self.output = [x.split(u'\n')[0] for x in out_data]
            except IOError:
                self.output = []
        elif isinstance(output_src, list):
            self.output = output_src
        else:
            raise ValueError(u"Invalid input type for class parser")

    #-----------------------------------------------------------------------
    def flagError(self, msg):
        """ Indicate an error was encountered.
        """
        self.err = True
        self.err_msg = str(msg)

    #-----------------------------------------------------------------------
    def getLineIndices(self, key, str_list=None):
        """
        Get the indices of a list of strings in which a string contains a key.

        Args:
            key (str): The key to look for.
            str_list (list of str): The list to search (default None --> self.output).

        Returns:
            list of int
        """
        if not str_list: str_list = self.output
        ind = [ i for i, x in enumerate(str_list) if x.find(key) > -1 ]
        return ind

    #-----------------------------------------------------------------------
    def notEmpty(self, key):
        """
        Check if a string is empty (or only spaces).

        Args:
            key (str): The string to check.

        Returns:
            bool
        """
        return key != u'' and not key.isspace()

    #-----------------------------------------------------------------------
    def breakLine(self, line):
        """
        Split a string on spaces and remove empty elements.

        Args:
            line (str): The line to format.

        Returns:
            list of str
        """
        if line.find(u'\n') > -1:
            currentLine  = line.split(u'\n')[0]
        else:
            currentLine = line
        brokenLine = currentLine.split()
        populated_pieces = list(filter(self.notEmpty, brokenLine))
        return populated_pieces

    #-----------------------------------------------------------------------
    def getAllNumbers(self, input_line):
        """
        Get a list of every number in a string using regular expressions.

        Args:
            input_line (str): The string to search.

        Returns:
            list of float
        """
        return [float(num) for num in re.findall(parser.reg_ex, input_line)]

    #-----------------------------------------------------------------------
    def getNumbers(self, input_line, position_in_line, exceptions=False, err_val=None):
        """
        Get a float from a string using regular expressions.

        Args:
            input_line (str): The string to search.
            position_in_line (int): The index of the list generated by splitting the
                input_line on spaces.
            exceptions (bool): If True, catch exceptions, otherwise let them be raised (default False).
            err_val (user defined): Value to return if parsing error was encountered.

        Returns:
            float, err_val
        """
        try:
            formatted_line = self.breakLine(input_line)
            numbers = re.findall(parser.reg_ex, formatted_line[position_in_line])
            if len(numbers) > 1 :
                raise IOError(u"Multiple values found when expected one.")
            number = float(numbers[0])
        except Exception as ex:
            if not exceptions:
                raise ex
            err_msg = u"In parser.getNumbers, an exception of type {0} occured. "+\
                      u"Arguments:\n{1!r}"
            message = err_msg.format(type(ex).__name__, ex.args)
            number  = self.float_err
            if err_val is not None: number = err_val
            print(message)
        return number

