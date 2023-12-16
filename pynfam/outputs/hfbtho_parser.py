# -------- Backwards Compatibility ---------
from __future__ import print_function
from __future__ import division
# -------------- Utilities -----------------
import numpy as np
# ----------- Relative Imports -------------
from .parser import parser
# ------------------------------------------

__version__ = u'2.0.0'
__date__    = u'2019-07-26'

#===============================================================================#
#                             CLASS hfbthoParser                                #
#===============================================================================#
class hfbthoParser(parser):
    """
    An hfbthoParser is a subclass of the parser class that contains an
    hfbtho output and the methods to parse it for values.

    Note that an hfbtho output for internal blocking contains multiple single
    run outputs. In this case parsed values from the blocking candidate with the
    lowest binding energy are returned.

    Args:
        hfbtho_output (str, list): A str is assumed to be path to output file, which
            is then read in and converted to a list of strings. Newline symbols
            are stripped.

    Attributes:
        keys (dict): A collection of strings to search for in the output.
        blocking_loop (list): If not empty, indicates blocking is active.
        energy (float): tEnergy (minimum tEnergy if internal blocking is on).
        imin (int): Index of internal blocking section with minimum tEnergy.
    """

    def __init__(self, hfbtho_output):
        parser.__init__(self, hfbtho_output)
        # Keys to look for
        self.keys = {u'reg_stage':u'### REGULAR STAGE',
                u'init_stage'    :u'### INITIAL STAGE',
                u'blocking_loop' :u'Initialization for the even-even core',
                u'conv_yes'      :u'iteration converged',
                u'conv_no'       :u'iterations limit interrupt',
                u'blocking_cands':u'Blocking candidates are',
                u'blocked_n'     :u'neutron Blocking:',
                u'blocked_p'     :u'proton  Blocking:',
                u'tenergy'       :u'tEnergy',
                u'def'           :u'deformation beta2',
                u'lam'           :u'lambda (ala)',
                u'time'          :u'Total CPU time',
                u'errors'        :[u'ERRORS IN HFBTHO_SOLVER',u'Segmentation fault'],
                u'qpn_start'     :u'#quasiparticle energies neutrons',
                u'qpp_start'     :u'#quasiparticle energies protons',
                u'qpn_end'       :u'#canonical s.p. energies neutrons',
                u'qpp_end'       :u'#canonical s.p. energies protons',
                u'sp_end'        :u'Sum canonical',
                u'cand_not_found' :u'No blocking candidate found'}
        # Internal blocking traits
        self.blocking_loop = self.getLineIndices(self.keys[u'blocking_loop'])
        self.energy = None
        self.imin = 0

    #-----------------------------------------------------------------------
    def findErrors(self):
        """ Look for internal HFBTHO errors. Set error attribute accordingly.
        """
        err_lines = []
        for err in self.keys[u'errors']:
            err_lines += self.getLineIndices(err)
        if err_lines:
            self.flagError(u"HFBTHO output indicates an internal error.")

        # Treat candidate not found explicitly
        if self.getLineIndices(self.keys[u'cand_not_found']):
            self.flagError(self.keys[u'cand_not_found'])

    #-----------------------------------------------------------------------
    def getBlockingCands(self):
        """
        Extract a list of blocking candidate quantum numbers.

        Returns:
            tuple of list

        Raises:
            IOError
        """
        try:
            inds = self.getLineIndices(self.keys[u'blocking_cands'])
            if not inds:
                raise IOError(u"No blocking candidates found in HFBTHO output.")
            cand_list = ([],[])
            # Only take candidates from the 1st soln in output (i.e. even core)
            for np_ind, line_ind in enumerate(inds[0:2]):
                counter = 2
                line = self.output[line_ind+counter]
                while self.notEmpty(line):
                    blocking_cand = 5*[0]
                    # Extract parity
                    halved = line.split(u'[')
                    sign_str = halved[0].split()[-1][-1]
                    if sign_str == u'-':
                        blocking_cand[1] = -1
                    elif sign_str == u'+':
                        blocking_cand[1] = 1
                    # Extract other quantum numbers
                    quant = self.getAllNumbers(line)[-4:]
                    blocking_cand[0] = int(quant[0]) # 2*abs(Omega)
                    blocking_cand[2] = int(quant[1]) # N
                    blocking_cand[3] = int(quant[2]) # nz
                    blocking_cand[4] = int(quant[3]) # nl
                    # Store candidate and update for next line
                    cand_list[np_ind].append(blocking_cand)
                    counter += 1
                    line = self.output[line_ind+counter]
        except Exception as e:
            cand_list = []
            self.flagError(e)

        return cand_list

    #-----------------------------------------------------------------------
    def getConv(self):
        """
        Determine whether all, some, or no solutions converged.

        Returns:
            str

        Raises:
            IOError
        """
        try:
            # Deal with kickoff mode...
            start_reg  = self.getLineIndices(self.keys[u'reg_stage'])
            start_init = self.getLineIndices(self.keys[u'init_stage'])
            if not start_init:
                # If kickoff mode off, just search
                conv_yes = self.getLineIndices(self.keys[u'conv_yes'])
                conv_no  = self.getLineIndices(self.keys[u'conv_no'])
            else:
                # If kickoff mode on, ignore lines for initial stage. We do this
                # by searching only blocks of data between the first occurance of
                # regular stage and the next occurence of initial stage.
                conv_yes, conv_no = [], []
                for i, start in enumerate(start_reg):
                    if i+1 == len(start_reg):
                        end = None
                    else:
                        end = start_init[i+1] # recall slice goes to end-1
                    sub_output = self.output[start:end]
                    conv_yes += self.getLineIndices(self.keys[u'conv_yes'], sub_output)
                    conv_no  += self.getLineIndices(self.keys[u'conv_no'], sub_output)

            # If internal blocking off, there should only be 1 occurence of conv
            if not self.blocking_loop:
                if len(conv_yes) + len(conv_no) > 1:
                    raise IOError(u"Unexpected convergence lines found in HFBTHO output.")
            if conv_yes:
                conv = u'Yes'
                # If internal blocking on, make sure ALL candidates converged
                if self.blocking_loop:
                    if conv_no:
                        conv = u'ATTN!'
            elif conv_no:
                conv = u'No'
            else:
                raise IOError(u"No convergence lines found in HFBTHO output.")
        except Exception as e:
            conv = self.str_err
            self.flagError(e)

        return conv

    #-----------------------------------------------------------------------
    def getEnergy(self):
        """
        Get the total binding energy.

        Returns:
            float

        Raises:
            IOError
        """
        try:
            energy_lines = self.getLineIndices(self.keys[u'tenergy'])
            energies = [self.getNumbers(self.output[line], -1) for line in energy_lines]
            if not energies:
                raise IOError(u"No total energy found in HFBTHO output.")
            self.imin = energies.index(min(energies))
            self.energy = min(energies)
        except Exception as e:
            self.imin = 0
            self.energy = self.float_err
            self.flagError(e)

        return self.energy

    #-----------------------------------------------------------------------
    def getTime(self):
        """
        Get the total run time.

        Returns:
            float
        """
        try:
            time_lines = self.getLineIndices(self.keys[u'time'])
            times = [self.getNumbers(self.output[line], -2) for line in time_lines]
            time = sum(times)
        except Exception as e:
            time = self.float_err
            self.flagError(e)

        return time

    #-----------------------------------------------------------------------
    def getLambda(self):
        """
        Get the proton and neutron fermi energies.

        Returns:
            tuple of float

        Raises:
            IOError
        """
        try:
            # Internal blocking: We have to run getEnergy to get most bound index
            if self.blocking_loop and self.energy is None:
                self.getEnergy()
            lambda_lines = self.getLineIndices(self.keys[u'lam'])
            if not lambda_lines:
                raise IOError(u"No fermi energy found in HFBTHO output.")
            # (neutrons, protons)
            lambda_vals = [(self.getNumbers(self.output[line], -2),
                            self.getNumbers(self.output[line], -1)) for line in lambda_lines]
        except Exception as e:
            lambda_vals = [2*(self.float_err,)]
            self.flagError(e)

        return lambda_vals[self.imin]

    #-----------------------------------------------------------------------
    def getBeta2(self):
        """
        Get the deformations beta2 of the final solution.

        Returns:
            tuple: (proton, neutron, total).

        Raises:
            IOError
        """
        try:
            # Internal blocking: We have to run getEnergy to get most bound index
            if self.blocking_loop and self.energy is None:
                self.getEnergy()
            beta2_lines = self.getLineIndices(self.keys[u'def'])
            if not beta2_lines:
                raise IOError(u"No deformation beta2 found in HFBTHO output.")
            # (neutrons, protons, total)
            beta2_vals = [(self.getNumbers(self.output[line], -3),
                           self.getNumbers(self.output[line], -2),
                           self.getNumbers(self.output[line], -1)) for line in beta2_lines]
        except Exception as e:
            beta2_vals = [3*(self.float_err,)]
            self.flagError(e)

        return beta2_vals[self.imin]

    #-----------------------------------------------------------------------
    def getSPData(self):
        """
        Parse for single particle data block.

        Returns:
            dict: {'p':prot_data_dict, 'n':neut_data_dict}.

        Raises:
            IOError
        """
        # Initialize p-n dict
        pn_data = {u'p': {}, u'n': {}}

        # Internal blocking: We have to run getEnergy to get most bound index
        if self.blocking_loop and self.energy is None: self.getEnergy()

        # Search for relevant lines
        pn_data[u'n'][u'start'] = self.getLineIndices(self.keys[u'qpn_end'])
        pn_data[u'p'][u'start'] = self.getLineIndices(self.keys[u'qpp_end'])
        pn_data[u'n'][u'end'] = self.getLineIndices(self.keys[u'sp_end'])
        pn_data[u'p'][u'end'] = self.getLineIndices(self.keys[u'sp_end'])

        try:
            if not (pn_data[u'n'][u'start'] and pn_data[u'n'][u'end']\
                and pn_data[u'p'][u'start'] and pn_data[u'p'][u'end']):
                raise IOError(u"Could not find SP data block in HFBTHO output.")

            # Extract the raw data per line into a list of dictionaries
            for t, pn in list(pn_data.items()):
                end = 0
                if t==u'p': end = 1
                qp_data  = self.output[pn[u'start'][self.imin]+10: pn[u'end'][2*self.imin+end]-1]
                qp_dicts = \
                [{u'ind'    : int(self.getNumbers(line,0)),
                  u'ce'     : self.getNumbers(line,2),
                  u'occ'    : self.getNumbers(line,3),
                  u'ovlp'   : self.getNumbers(line,6),
                  u'omega'  : int(self.getAllNumbers(line)[-4:][0]),
                  u'N'      : int(self.getAllNumbers(line)[-4:][-3]),
                  u'nz'     : int(self.getAllNumbers(line)[-4:][-2]),
                  u'nl'     : int(self.getAllNumbers(line)[-4:][-1]),
                  u'par'    : int(((line.split(u'[')[0]).split()[-1][-1])+u'1'),
                  u'ns'     : int(self.getAllNumbers(line)[-4:][0]) - 2*int(self.getAllNumbers(line)[-4:][-1])
                  } for line in qp_data]
                pn[u'qp_dicts'] = qp_dicts
        except Exception as e:
            pn_data = None
            self.flagError(e)

        return pn_data

    #-----------------------------------------------------------------------
    def getQPData(self):
        """
        Parse for quasi-particle data block.

        Returns:
            dict: {'p':prot_data_dict, 'n':neut_data_dict}.

        Raises:
            IOError
        """
        # Initialize p-n dict
        pn_data = {u'p': {}, u'n': {}}

        # Internal blocking: We have to run getEnergy to get most bound index
        if self.blocking_loop and self.energy is None: self.getEnergy()

        # Search for relevant lines
        pn_data[u'n'][u'start'] = self.getLineIndices(self.keys[u'qpn_start'])
        pn_data[u'p'][u'start'] = self.getLineIndices(self.keys[u'qpp_start'])
        pn_data[u'n'][u'end'] = self.getLineIndices(self.keys[u'qpn_end'])
        pn_data[u'p'][u'end'] = self.getLineIndices(self.keys[u'qpp_end'])
        times = self.getLineIndices(self.keys[u'time'])
        start_reg = self.getLineIndices(self.keys[u'reg_stage'])

        try:
            if not (times and start_reg and\
                    pn_data[u'n'][u'start'] and pn_data[u'n'][u'end'] and\
                    pn_data[u'p'][u'start'] and pn_data[u'p'][u'end']):
                raise IOError(u"Could not find QP data block in HFBTHO output.")

            # Extract the raw data per line into a list of dictionaries
            for t, pn in list(pn_data.items()):
                qp_data  = self.output[pn[u'start'][self.imin]+8: pn[u'end'][self.imin]-5]
                qp_dicts = \
                [{u'ind'    : int(self.getNumbers(line,0)),
                  u'block'  : int(self.getNumbers(line,1)),
                  u'eqp'    : self.getNumbers(line,2),
                  u'occ'    : self.getNumbers(line,6),
                  u'overl'   : self.getNumbers(line,8),
                  u'omega'  : int(self.getAllNumbers(line)[-4:][0]),
                  u'N'      : int(self.getAllNumbers(line)[-4:][-3]),
                  u'nz'     : int(self.getAllNumbers(line)[-4:][-2]),
                  u'nl'     : int(self.getAllNumbers(line)[-4:][-1]),
                  u'par'    : int(((line.split(u'[')[0]).split()[-1][-1])+u'1'),
                  u'ns'     : int(self.getAllNumbers(line)[-4:][0]) - 2*int(self.getAllNumbers(line)[-4:][-1])
                  } for line in qp_data]

                # Get blocking candidate info from below iterations printout
                blo_dict = None
                if t == u'p':
                    blo_on = self.getLineIndices(self.keys[u'blocked_p'])
                else:
                    blo_on = self.getLineIndices(self.keys[u'blocked_n'])
                if blo_on:
                    # Careful to skip kickoff mode time
                    itime = times[np.argmax(np.array(times) > start_reg[self.imin])]
                    b_line_ind = blo_on[np.argmax(np.array(blo_on) > itime)]
                    blo_line   = self.output[b_line_ind]
                    blo_dict = {u'eqp'  : self.getNumbers(blo_line,7),
                                u'block': int(self.getNumbers(blo_line,3)),
                                u'state': int(self.getNumbers(blo_line,5)) - 1,
                                u'omega': int(self.getAllNumbers(blo_line)[-4:][0]),
                                u'N'    : int(self.getAllNumbers(blo_line)[-4:][-3]),
                                u'nz'   : int(self.getAllNumbers(blo_line)[-4:][-2]),
                                u'nl'   : int(self.getAllNumbers(blo_line)[-4:][-1]),
                                u'par'  : int(((blo_line.split(u'[')[0]).split()[-1][-1])+u'1')
                               }
                    blo_dict[u'ns'] = blo_dict[u'omega'] - 2*blo_dict[u'nl']
                pn[u'qp_dicts'] = qp_dicts
                pn[u'blo_qp']   = blo_dict
        except Exception as e:
            pn_data = None
            self.flagError(e)

        return pn_data
