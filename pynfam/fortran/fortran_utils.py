# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
from builtins   import object
# -------------- Utilities -----------------
from shutil import copy2
import os
import stat
import sys
import copy
import f90nml
"""
Use subprocess32 in POSIX system for Python 2 as subprocess in Python 2 is not OFED-fork-safe.
In Python 3 subprocess is replaced by subprocess32 for POSIX system so no need to bother.
Ref: https://bitbucket.org/fenics-project/instant/pull-requests/13/try-to-import-subprocess32-fall-back-to/diff
Ref: https://fenics.readthedocs.io/projects/instant/en/latest/installation.html
Ref: https://www.open-mpi.org/faq/?category=openfabrics#ofa-fork
"""
if os.name == u"posix" and sys.version_info[0] < 3:
    try:
        import subprocess32 as subprocess
    except ImportError:
        print(u"Warning: subprocess32 cannot be imported.")
        print(u"Warning: subprocess32 is strongly recommended as subprocess is not OFED-fork-safe!")
        import subprocess
else:
    import subprocess
# ----------- Relative Imports -------------
from ..outputs.pynfam_paths import pynfamPaths
# ------------------------------------------
from glob import glob

__version__ = u'2.0.0'
__date__    = u'2019-07-26'

#===============================================================================#
#                             CLASS nmlInputs                                   #
#===============================================================================#
class nmlInputs(object):
    """
    An instance of the class nmlInputs contains a deep copy of a fortran namelist
    (as an ordered dict) and the methods to read, write and adjust parameter
    values using f90nml.

    This is a base class for the fortProcess class.

    Args:
        fname (str): The namelist filename.
        nml (OrderedDict): The namelist (default None).

    Attributes:
        fname_nml (str): The namelist filename.
        nml (OrderedDict): The namelist (default None).
    """

    def __init__(self, fname, nml=None):
        self.fname_nml = fname
        if nml is not None:
            self.nml = copy.deepcopy(nml) # Ordered Dict
        else:
            self.nml = None

    #-----------------------------------------------------------------------
    def readNml(self, source):
        """
        Populate nml attribute with an f90nml object (an ordered dict) read
        from an existing namelist text file.

        Args:
            source (str): Path to namelist file.
        """

        f2open = os.path.join(source, self.fname_nml)
        with open(f2open, u'r') as nml_file:
            nml_obj = f90nml.read(nml_file)
        self.nml = nml_obj

    #-----------------------------------------------------------------------
    def writeNml(self, dest=u'./'):
        """
        Convert ordered-dict to f90nml object and write to file.

        Args:
            dest (str): Path where namelist file is written (default './').
        """

        dest = os.path.join(dest, self.fname_nml)

        HFBTHO_nml = f90nml.Namelist(self.nml)
        with open(dest, u'w') as nml_file:
            HFBTHO_nml.write(nml_file)

    #-----------------------------------------------------------------------
    def adjustNml(self, inputs):
        """
        Adjust the values in the namelist.

        Accepts a 2 level dictionary mapping to the namelist that only contains
        the parameters to be changed. If the namelist value is a list, provide
        None for list elements which should remain unchanged.

        Args:
            inputs (dict): Two level dictionary mapping to namelist entries.
        """

        for h1, d1 in list(inputs.items()):
            for h2, param in list(d1.items()):
                # if we need to adjust a list, only adjust active elements
                if isinstance(param, list):
                    self.nml[h1][h2] = [p if p is not None else self.nml[h1][h2][i]\
                            for i,p in enumerate(param)]
                else:
                    # Allow passing None or '' to get fortran code's default.
                    # f90nml turns None into blank line.
                    if param is None:
                        self.nml[h1][h2] = None
                    elif param == u'':
                        self.nml[h1][h2] = None
                    else:
                        self.nml[h1][h2] = param

    #-----------------------------------------------------------------------
    def setNmlParam(self, param_dict):
        """
        Convenience method to adjust namelist parameter with a 1 level dict of
        the form {'param': value} (opposed to {'header':{'param':value}}).

        Args:
            param_dict (dict): One level dict with value to change.
        """

        header_dict = {h:[hh for hh in self.nml[h]] for h in self.nml}

        # Check we don't have any extraneous params
        all_params = [p for h in header_dict for p in header_dict[h]]
        rogues     = [p for p in list(param_dict.keys()) if p not in all_params]
        if rogues:
            raise ValueError(u"Invalid namelist input detected: {:}".format(rogues))

        # Group params by header and construct 2 tier dict
        in_dict = {}
        for h in header_dict:
            h_dict = {}
            for p in param_dict:
                if p in header_dict[h]: h_dict[p] = param_dict[p]
            in_dict[h] = h_dict

        # Adjust the namelist accordingly
        self.adjustNml(in_dict)


#===============================================================================#
#                             CLASS fortProcess                                 #
#===============================================================================#
class fortProcess(nmlInputs):
    """
    A fortProcess is a subclass of the nmlInputs class, with additional attributes
    and method for running a fortran executable in a self-contained directory.

    Args:
        paths (pynfamPaths): The paths for the pynfam calculation.
        exe (str): Fortran executable file name.
        fname_nml (str): The namelist filename.
        nml (OrderedDict): The namelist (default None).
        label (str): Unique name for run directory (default None).

    Attributes:
        paths (pynfamPaths): The paths for the pynfam calculation.
        exe (str): Fortran executable file name.
        label (str): Unique name for run directory.
        output (list of str): Lines of the fortran program's output (excluding
            newline charcters).
    """

    def __init__(self, paths, exe, fname_nml, nml=None, label=None):
        if not isinstance(paths, pynfamPaths):
            raise TypeError(u"Paths must be a pynfamPaths object")
        self.paths = copy.deepcopy(paths)
        self.exe   = exe
        self.label = label
        self.tmp_inputs = []
        self.pmt_inputs = []
        self.paths.exe = self.find_exe(self.paths.exe, self.exe)
        if self.paths.exe is None:
            raise IOError(u"Fortran executable "+exe+" not found")
        nmlInputs.__init__(self, fname_nml, nml)

    #-----------------------------------------------------------------------
    @staticmethod
    def find_exe(exe_root, exe):
        """ 
        Find the Fortran executable (filename given by exe) in directory 
        exe_root and its subdirectories. Assume there is only one file 
        with the given filename in directory exe_root and its subdirectories.
        Execute permission will be granted if it does not have it.

        Returns:
            None, str: Path of the directory that contains the Fortran executable,
                       or None if such executable is not found.
        """
      
        if not os.path.isdir(exe_root):
            return None
        file_path = None
        for dirpath, _, filenames in os.walk(exe_root, topdown=True):
            if exe in filenames:
                file_path = os.path.join(dirpath, exe)
                break
        if file_path is None:
            return None
        if not os.access(file_path, os.X_OK):
            os.chmod(file_path, os.stat(file_path).st_mode | stat.S_IEXEC)
        return dirpath

    @property
    def rundir(self):
        """ rundir (str): The full path to run directory.
        """
        return os.path.join(self.paths.out,self.label)

    #-----------------------------------------------------------------------
    def runExe(self, stdout_bool=False, stdout_return=False, statement=None, debug=0):
        """
        Launch a fortran executable from rundir and capture stdout and stderr.
        Pipe captured stdout to sys stdout in real time if desired.

        Args:
            stdout_bool (bool): Print to stdout in real time (default False).
            stdout_return (bool): Return stdout (default False).
            statement (str or list): Statement for Popen with shell=False
                (default None --> '/path2exe/exe').
            debug (int): If !=0, don't run the exe (default 0).

        Returns:
            list      : Stdout as list of strings; empty if stdout_return is False.
            None, str : Stderr, or None if stderr is empty.
        """

        # Setup directory to run the program in
        if not self.label:
            raise RuntimeError(u"Cannot execute run_exe until run label has been set")
        if not os.path.exists(self.rundir):
            os.makedirs(self.rundir)

        # Populate rundir with the necessary inputs files
        self.writeNml(self.rundir) # namelist
        pmt_inputs_dest = [copy2(other, self.rundir) for other in self.pmt_inputs] # permanent
        tmp_inputs_dest = [copy2(other, self.rundir) for other in self.tmp_inputs] # temporary

        # Default execution statement is to run exe from current directory
        if statement is None:
            statement = [os.path.join(self.paths.exe,self.exe)]

        # Run the program
        os.chdir(self.rundir)

        output, err = [], u''

        if debug == 0:
            if stdout_bool:
                stdout_handle = subprocess.PIPE
                stderr_handle = subprocess.PIPE
            else:
                stderr_fname = 'stderr.tmp'
                stderr_handle = open(stderr_fname, 'w+')
                if stdout_return:
                    # temp files for stdout and stderr
                    stdout_fname = 'stdout.tmp'
                    stdout_handle = open(stdout_fname, 'w+')
                else:
                    stdout_handle = subprocess.DEVNULL
            # Popen docs: 'universal_newlines' is equivalent to 'text' in versions >~ 3,
            # and is provided for backwards compatibility.
            proc = subprocess.Popen(statement, shell=False, universal_newlines=True,
                    stdin=subprocess.DEVNULL, stdout=stdout_handle, stderr=stderr_handle)
            if stdout_bool:
                # Print in real time (proc.stdout is a file object)
                while True:
                    so = proc.stdout.readline()
                    if so == u'' and proc.poll() is not None:
                        break
                    if so:
                        sys.stdout.write(so)
                        sys.stdout.flush()
                        if stdout_return:
                            output.append(so.split(u'\n')[0])
                err = proc.stderr.readline()
            else:
                # Wait until proc is finished and get stdout
                # output_data, err = proc.communicate()
                # output = [l for l in output_data.split(u'\n')]
                proc.wait()
                if stdout_return:
                    stdout_handle.seek(0)
                    output = stdout_handle.readlines()
                    stdout_handle.close()
                stderr_handle.seek(0)
                err = stderr_handle.readline()
                stderr_handle.close()

        os.chdir(self.paths.cwd)

        # Obtain stderr (a string or None)
        if err == u'' or err.isspace():
            self.delete_tmps(tmp_inputs_dest) # Delete temporary files
            err = None

        # This is wrapped by sub-classes which only return err but use output
        return output, err

    def delete_tmps(self, tmp_inputs_dest):
        """
        Delete temporary input and output files.

        Args:
            tmp_inputs_dest: List of temporary input files.
        """
        for x in tmp_inputs_dest + glob(os.path.join(self.rundir, u'*.tmp')):
            if os.path.isfile(x):
                os.remove(x)