# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
from builtins   import object
from builtins   import str
# -------------- Utilities -----------------
from shutil import copy2
import os
import errno
#from pathlib2 import Path
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-07-26'

#===============================================================================#
#                             CLASS pynfamPaths                                 #
#===============================================================================#
class pynfamPaths(object):
    """
    An instance of the class pynfamPaths contains absolute paths to directories
    comprising the output directory tree for a single pynfam calc. It has methods
    to create the directory tree, parse if for labels and copy and delete files.
    The tree is structed as follows, where parenthesis indicate files, and *
    indicates directories not made by this class::

        outputs
          |- meta
          |    |- (copied_inputs_and_slurm_files)
          |
          |- calc_label
               |- hfb_soln
               |    |- (copied_ground_state_output_files)
               |    |- hfb_meta
               |         |- (logfile_for_all_hfb_runs)
               |         |- (tarfile_of_all_hfb_run_dirs)
               |              |- hfb_label_even*
               |              |    |- (raw_output_files)
               |              |- hfb_label_odd*
               |                   |- core*
               |                       |- (raw_output_files)
               |                   |- candidate_label*
               |                       |- (raw_output_files)
               |
               |- fam_soln
               |    |- (strength_output_text_and_binary_files_per_operator)
               |    |- fam_meta
               |         |- (representative_logfile_for_pnfam_runs)
               |         |- (tarfiles_of_all_pnfam_runs_per_operator)
               |              |- fam_point_label*
               |                   |- (raw_output_files)
               |
               |- beta_soln
                    |- (beta_decay_output_files)
                    |- beta_meta
                         |- (non_essential_output_files)

    Args:
        ls_dir (str): name of the top level output directory.
        exe_dir (str): name of the directory containing executables.
        scr_dir (str): name of the scratch working directory.

    Attributes:
        cwd (str): absolute path to the current working directory.
        exe (str): absolute path to the executables.
        scr (str): absolute path to the scratch working directory.
        calclabel (str): unique label for pynfam directory under the top level directory..
        _lsdir (str): name of the top level output directory.
    """

    subdir_meta  = u'meta'
    """ subdir_meta (str): subdirectory name for easy maintentance.
    """
    subdir_hfb   = u'hfb_soln'
    """ subdir_hfb (str): subdirectory name for easy maintentance.
    """
    subdir_hfb_m = u'hfb_meta'
    """ subdir_hfb_m (str): subdirectory name for easy maintentance.
    """
    subdir_fam   = u'fam_soln'
    """ subdir_fam (str): subdirectory name for easy maintentance.
    """
    subdir_fam_m = u'fam_meta'
    """ subdir_fam_m (str): subdirectory name for easy maintentance.
    """
    subdir_beta  = u'beta_soln'
    """ subdir_beta (str): subdirectory name for easy maintentance.
    """
    subdir_beta_m= u'beta_meta'
    """ subdir_beta_m (str): subdirectory name for easy maintentance.
    """
    subdir_hfbcore = u'core'
    """ subdir_hfbcore (str): subdirectory name for easy maintentance.
    """
    ignore = [u'meta',u'soln']
    """ ignore (list of str): Keys to ignore in directory searchs.
    """

    def __init__(self, label, ls_dir=u'./', exe_dir=u'./', scr_dir=u'./'):
        self.cwd = os.getcwd()
        self.exe = os.path.abspath(exe_dir)
        self.scr = os.path.abspath(scr_dir)
        self._lsdir = ls_dir
        self.calclabel = str(label)
        if label is None:
            self.calclabel = label

    @property
    def out(self):
        """ out (str): absolute path to 'scr/_lsdir'.
        """
        return os.path.join(self.scr, self._lsdir)
    @property
    def meta(self):
        """ meta (str): absolute path to 'out/subdir_meta'.
        """
        return os.path.join(self.out, pynfamPaths.subdir_meta)
    @property
    def calc(self):
        """ calc (str): absolute path to 'out/calclabel'.
        """
        return os.path.join(self.out, self.calclabel)
    @property
    def hfb(self):
        """ hfb (str): absolute path to 'calc/subdir_hfb'. Requires calclabel!=None.
        """
        return os.path.join(self.calc, pynfamPaths.subdir_hfb)
    @property
    def fam(self):
        """ fam (str): absolute path to 'calc/subdir_fam'. Requires calclabel!=None.
        """
        return os.path.join(self.calc, pynfamPaths.subdir_fam)
    @property
    def beta(self):
        """ beta (str): absolute path to 'calc/subdir_beta'. Requires calclabel!=None.
        """
        return os.path.join(self.calc, pynfamPaths.subdir_beta)
    @property
    def hfb_m(self):
        """ hfb_m (str): absolute path to 'hfb/subdir_hfb_m'. Requires calclabel!=None.
        """
        return os.path.join(self.hfb, pynfamPaths.subdir_hfb_m)
    @property
    def fam_m(self):
        """ fam_m (str): absolute path to 'fam/subdir_fam_m'. Requires calclabel!=None.
        """
        return os.path.join(self.fam, pynfamPaths.subdir_fam_m)
    @property
    def beta_m(self):
        """ beta_m (str): absolute path to 'beta/subdir_beta_m'. Requires calclabel!=None.
        """
        return os.path.join(self.beta, pynfamPaths.subdir_beta_m)

    #-----------------------------------------------------------------------
    def bn(self, full_path):
        """ Alias to get basename of an absolute path.
        """
        return os.path.basename(os.path.normpath(full_path))

    #-----------------------------------------------------------------------
    def rp(self, full_path):
        """ Alias to get the relative path from out.
        """
        return os.path.relpath(full_path, self.out)

    #-----------------------------------------------------------------------
    def addTag(self, tag):
        """ Add a tag to the end of the out directory.
        """
        self._lsdir = self.bn(self.out) + str(u'_'+tag)

    #-----------------------------------------------------------------------
    def mkdirs(self, newpath=None):
        """
        MPI compatible directory creation. Sometimes even when only rank = 0
        calls os.mkdir() it still raises errors on the other processes and
        causes the program to exit. This avoids such a scenario.
        """
        # If no argument, make as many directories in the tree as available
        if newpath is None:
            if not self.calclabel:
                paths2mk = [self.out, self.meta]
            else:
                paths2mk = [self.out, self.meta, self.calc,
                            self.hfb, self.hfb_m,
                            self.fam, self.fam_m,
                            self.beta,self.beta_m]
        # Otherwise assume the argument is a full path to make
        elif isinstance(newpath, list):
            paths2mk = newpath
        else:
            paths2mk = [newpath]

        # Try making the directory, skipping those raised if already exists
        # (e.g. by other tasks)
        while True:
            try:
                for path in paths2mk:
                    if not os.path.exists(path):
                        os.mkdir(path)
                break
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise
                pass

    #-----------------------------------------------------------------------
    def getLabels(self, top=None, tier2=False, prefix=None):
        """
        Get sorted directory names, excluding those with 'ignore' keywords.

        * Optional 1 tier or 2 tier label

            - e.g. tier2=True, top=out: 000000
            - e.g. tier2=True, top=fam_m: F-K0/000000 or 0000/core

        * Optional return of fullpath or relative path from specified prefix. If
          prefix=True, return fullpath, if prefix=string, return relative path from string.

            - e.g. prefix=True, top=hfb_m: /path2top/top/000000/hfb_soln/hfb_meta/0000/core
            - e.g. prefix=out,  top=fam_m: 000000/fam_soln/fam_meta/F-K0/000000
        """
        if top is None: top = self.out

        alllabels = []
        dir_contents = os.listdir(top)
        for top_label in dir_contents:
            fullname = os.path.join(top, top_label)
            ignoreit = any(ign in top_label for ign in pynfamPaths.ignore)
            if os.path.isdir(fullname) and not ignoreit:
                labels2 = []
                label1_too = []
                if tier2:
                    dir_contents2 = os.listdir(fullname)
                    for next_label in dir_contents2:
                        fullname2 = os.path.join(fullname, next_label)
                        ignoreit2 = any(ign in next_label for ign in pynfamPaths.ignore)
                        if os.path.isdir(fullname2) and not ignoreit2:
                            labels2.append(next_label)
                        # If non-directory contents in toplabel, store it as a unique label
                        if not os.path.isdir(fullname2) and not ignoreit2:
                            label1_too.append(next_label)
                if labels2:
                    alllabels += [os.path.join(top_label, l2) for l2 in labels2]
                if not labels2 or label1_too:
                    alllabels += [top_label]

        if not prefix: # None or False or Empty is fine
            finallabels = alllabels
        elif prefix is True: #Explicit check for boolean here
            finallabels = [os.path.join(top, l) for l in alllabels]
        elif isinstance(prefix, str):
            finallabels = [os.path.relpath(os.path.join(top,l), prefix) for l in alllabels]
        else:
            raise TypeError(u"Invalid argument type for kwarg prefix.")

        return sorted(finallabels)

    #-----------------------------------------------------------------------
    def copyAllFiles(self, src, dest):
        """ Copy all files (only) from src to dest directory.
        """
        src_files = os.listdir(src)
        for file_name in src_files:
            full_file_name = os.path.join(src, file_name)
            if (os.path.isfile(full_file_name)):
                copy2(full_file_name, dest)

    #-----------------------------------------------------------------------
    def deleteAll(self, src, dest):
        """ List all files (only) in src and delete them from dest.
        """
        src_files = os.listdir(src)
        file_list = [os.path.join(dest, f) for f in src_files]
        for f2delete in file_list:
            if os.path.isfile(f2delete):
                os.remove(f2delete)
