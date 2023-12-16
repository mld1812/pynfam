#!/usr/bin/env python
# Run a basic strength function test to check for regressions

import datetime
import os
import shlex
import subprocess as sp
import sys


class Operator(object):
    """Container for a PNFAM operator.
    
    Attributes:
        k           -- a list of K projections (e.g. [1, 0, -1])
        label       -- the PNFAM operator label (e.g. "GT")
        beta_minus  -- True if beta-, False if beta+
        fancy_label -- the PNFAM fancy label (e.g. "GT-")
    
    """
    
    def __init__(self, label, k):
        """Be smart and set class properties using arguments.
        
        Args:
            label -- a string label that must end in "+" or "-"
            k     -- an integer representing |max K| for the operator
        
        """
        self.k = sorted(range(-k, k+1, 1), key=lambda s: (abs(s), -s))
        label = label.strip().upper()
        if label[-1] not in '+-':
            raise ValeuError("Unknown beta+-")
        self.label      = label[:-1]
        self.beta_minus = True if label[-1] == '-' else False
    
    def __str__(self):
        return '{} (J={})'.format(self.fancy_label, max(self.k))
    
    @property
    def fancy_label(self):
        """Combine label and beta_minus attributes to form a fancy label."""
        return '{}{}'.format(self.label, '-' if self.beta_minus else '+')


## -----------------------------------------------------------------------------
## DEFINITIONS
## -----------------------------------------------------------------------------

STRENGTH_TOLERANCE = 1.0e-7

OPERATORS_MINUS = [Operator('F-', 0), Operator('GT-', 1), Operator('RS0-', 0),
   Operator('PS0-', 0), Operator('RS1-', 1), Operator('R-', 1),
   Operator('P-', 1), Operator('RS2-', 2)]

OPERATORS_PLUS = [Operator('F+', 0), Operator('GT+', 1), Operator('RS0+', 0),
   Operator('PS0+', 0), Operator('RS1+', 1), Operator('R+', 1),
   Operator('P+', 1), Operator('RS2+', 2)]

STRENGTH_DATA = dict()
STRENGTH_DATA['F-:0']    = [5.000000, 1.4873996454695729e-02]
STRENGTH_DATA['GT-:0']   = [5.000000, 3.4537324312349023e-02]
STRENGTH_DATA['GT-:1']   = [5.000000, 5.2172811191411782e-02]
STRENGTH_DATA['GT-:-1']  = [5.000000, 5.2172811191411754e-02]
STRENGTH_DATA['RS0-:0']  = [5.000000, 3.0797931541368919e-02, 3.0797931541368919e-02,-3.1825226921143354e-02]
STRENGTH_DATA['PS0-:0']  = [5.000000, 9.0861056404771979e-04,-3.1825213258481158e-02, 9.0861056404771979e-04]
STRENGTH_DATA['RS1-:0']  = [5.000000, 1.4315571166918362e-01,-5.0744589893595794e-02, 1.4315571166918362e-01, 1.0932110397100046e-02]
STRENGTH_DATA['RS1-:1']  = [5.000000, 7.1809808713095544e-02,-2.2598900430001501e-02, 7.1809808713095544e-02, 3.7949435463788240e-03]
STRENGTH_DATA['RS1-:-1'] = [5.000000, 7.1809808713094725e-02,-2.2598900430001401e-02, 7.1809808713094725e-02, 3.7949435463788326e-03]
STRENGTH_DATA['R-:0']    = [5.000000, 3.6304991959993656e-02, 3.6304991959993656e-02,-5.0744594328558178e-02, 5.7252972160243016e-03]
STRENGTH_DATA['R-:1']    = [5.000000, 3.3453873985566297e-02, 3.3453873985566297e-02,-2.2598906604974010e-02, 8.2228082547626657e-03]
STRENGTH_DATA['R-:-1']   = [5.000000, 3.3453873985566339e-02, 3.3453873985566339e-02,-2.2598906604973833e-02, 8.2228082547626518e-03]
STRENGTH_DATA['P-:0']    = [5.000000, 3.8693634115586142e-04, 5.7252924397958316e-03, 1.0932121281103515e-02, 3.8693634115586142e-04]
STRENGTH_DATA['P-:1']    = [5.000000, 7.5920413044835427e-05, 8.2228172825583833e-03, 3.7949389485016014e-03, 7.5920413044835427e-05]
STRENGTH_DATA['P-:-1']   = [5.000000, 7.5920413044827689e-05, 8.2228172825584214e-03, 3.7949389485016347e-03, 7.5920413044827689e-05]
STRENGTH_DATA['RS2-:0']  = [5.000000, 4.4292547757685985e-01]
STRENGTH_DATA['RS2-:1']  = [5.000000, 3.9937894107237387e-01]
STRENGTH_DATA['RS2-:-1'] = [5.000000, 3.9937894107237432e-01]
STRENGTH_DATA['RS2-:2']  = [5.000000, 3.4570766500850014e-01]
STRENGTH_DATA['RS2-:-2'] = [5.000000, 3.4570766500850070e-01]
STRENGTH_DATA['F+:0']    = [5.000000,-4.1798081988816084e-03]
STRENGTH_DATA['GT+:0']   = [5.000000, 5.1806214480011178e-03]
STRENGTH_DATA['GT+:1']   = [5.000000, 2.8714928484135115e-02]
STRENGTH_DATA['GT+:-1']  = [5.000000, 2.8714928484135083e-02]
STRENGTH_DATA['RS0+:0']  = [5.000000, 2.4504418360840186e-01, 2.4504418360840186e-01,-1.1896568654305310e-01]
STRENGTH_DATA['PS0+:0']  = [5.000000, 4.0911876920766047e-02,-1.1896567228962107e-01, 4.0911876920766047e-02]
STRENGTH_DATA['RS1+:0']  = [5.000000, 7.7677032755427922e+00, 2.6338487170300930e+00, 7.7677032755427922e+00, 1.2901834290713241e+00]
STRENGTH_DATA['RS1+:1']  = [5.000000, 5.6600302171384422e+00, 3.1109835500944927e+00, 5.6600302171384422e+00, 1.4758253968921105e+00]
STRENGTH_DATA['RS1+:-1'] = [5.000000, 5.6600302171384378e+00, 3.1109835500944905e+00, 5.6600302171384378e+00, 1.4758253968921120e+00]
STRENGTH_DATA['R+:0']    = [5.000000, 1.0772155226669078e+00, 1.0772155226669078e+00, 2.6338486444452647e+00, 5.3492972907198100e-01]
STRENGTH_DATA['R+:1']    = [5.000000, 2.0900687246861036e+00, 2.0900687246861036e+00, 3.1109836344335191e+00, 9.9626912032883397e-01]
STRENGTH_DATA['R+:-1']   = [5.000000, 2.0900687246861120e+00, 2.0900687246861120e+00, 3.1109836344335253e+00, 9.9626912032883341e-01]
STRENGTH_DATA['P+:0']    = [5.000000, 2.5925076607322328e-01, 5.3492973385007525e-01, 1.2901833820872137e+00, 2.5925076607322328e-01]
STRENGTH_DATA['P+:1']    = [5.000000, 4.6818809417872803e-01, 9.9626913714278031e-01, 1.4758255175892994e+00, 4.6818809417872803e-01]
STRENGTH_DATA['P+:-1']   = [5.000000, 4.6818809417872714e-01, 9.9626913714277920e-01, 1.4758255175892996e+00, 4.6818809417872714e-01]
STRENGTH_DATA['RS2+:0']  = [5.000000, 1.3410132374124666e+00]
STRENGTH_DATA['RS2+:1']  = [5.000000, 1.1853011433984471e+00]
STRENGTH_DATA['RS2+:-1'] = [5.000000, 1.1853011433984495e+00]
STRENGTH_DATA['RS2+:2']  = [5.000000, 3.5913491115688245e-01]
STRENGTH_DATA['RS2+:-2'] = [5.000000, 3.5913491115688384e-01]


def log(s):
    """Write string s to stderror with timestamp."""
    timestamp = datetime.datetime.strftime(datetime.datetime.now(), "(%H:%M)")
    sys.stderr.write('{}  {}\n'.format(timestamp, s.strip()))


def exec_pnfam(infile, outfile):
    """Run PNFAM and store the output.
    
    Args:
        infile  -- PNFAM namelist file
        outfile -- PNFAM output file (NOT the log file which goes to /dev/null)
    
    Returns the contents of outfile.
    
    """
    cmd = '../../pnfam_nompi.x "{}"'.format(infile)
    with open('/dev/null', 'w') as devnull:
        p = sp.Popen(shlex.split(cmd.format(infile)), stdout=devnull)
        retcode = p.wait()
    
    with open(outfile, 'r') as f:
        output = f.read()
    return output


def run_pnfam_str_calculations():
    """Run PNFAM calculations for all operators and compare to cached data.
    
    Returns number of tests run, number of tests with errors, and a
    dictionary of computed strength functions for testing rotational
    invariance.
    
    """
    with open('pnfam-str.tmpl', 'r') as f:
        tmpl = f.read()
    
    # Count tests & errors; save strength calculations for further testing
    tests  = 0
    errors = 0
    strengths = dict()
    
    for op in (OPERATORS_MINUS + OPERATORS_PLUS):
        for k in op.k:
            with open("test.in", 'w') as f:
                f.write(tmpl.format(label=op.fancy_label, k=k))
            
            # Run the code and get the strength --- on the last line of output
            output = exec_pnfam(infile='test.in', outfile='test.out')
            str_outs = output.strip().split('\n')[-1].split()
            str_head = output.strip().split('\n')[-2].split()[1:]
            str_vals = [v for v, h in zip(str_outs, str_head) if 'EQRPA' in h or 'Im' in h]
            strength = map(float, str_vals)
            
            key = '{}:{}'.format(op.fancy_label, k)
            strengths[key] = strength
            tests += 1
            
            # Equal shapes?
            if len(strength) != len(STRENGTH_DATA[key]):
                log("[{:<5s}{:>2d}]:  UNEQUAL SHAPES!".format(
                    op.fancy_label + ',', k))
                log(">> [{}] (calc)".format(", ".join(map(str, strength))))
                log(">> [{}] (expected)".format(
                    ", ".join(map(str, STRENGTH_DATA[key]))))
                log("")
                errors += 1
                continue
            
            # Diffs between calculation and cache
            diffs = map(abs, [strength[i]-STRENGTH_DATA[key][i]
                for i in xrange(len(strength))])
            
            try:
                if any([d > STRENGTH_TOLERANCE for d in diffs]):
                    raise
            except:
                errors += 1
                flag = "ERROR"
                continue
            else:
                flag = "ok"
            finally:
                log("[{:<5s}{:>2d}]:  [{}]  {}".format(op.fancy_label+',',
                    k, ", ".join(map(lambda s: '%.2e' % s, diffs)), flag))
    
    return tests, errors, strengths


def run_pnfam_str_rot_analysis(strength):
    """Check previously-computed strength (from run_pnfam_str_calculations)
    for invariance when K -> -K.
    
    Returns the number of tests and the number of tests which had errors.
    
    """
    tests  = 0
    errors = 0
    
    for op in (OPERATORS_MINUS + OPERATORS_PLUS):
        # No test when K = 0
        if max(op.k) > 0:
            for k in xrange(1, 1+max(op.k)):
                tests += 1
                
                key1 = '{}:{}'.format(op.fancy_label, k)
                key2 = '{}:{}'.format(op.fancy_label, -k)
                
                diffs = map(abs, [strength[key1][i]-strength[key2][i]
                    for i in xrange(len(strength[key1]))])
                
                if any([d > STRENGTH_TOLERANCE for d in diffs]):
                    errors += 1
                    flag = 'ERROR'
                else:
                    flag = 'ok'
                
                log("[{:<5s}{:>2d}]:  [{}]  {}".format(op.fancy_label+',', k,
                    ", ".join(map(lambda s: '%.2e' % s, diffs)), flag))
    
    return tests, errors


## -----------------------------------------------------------------------------
## MAIN PROGRAM
## -----------------------------------------------------------------------------

if __name__ == '__main__':
    # Make sure the HFB solution is here
    if not os.path.exists('solution.hfb'):
        log("ERROR: could not find file 'solution.hfb'")
        quit()
    
    # Run the strength calculations
    log("Running STR tests")
    log("")
    tests, errors, strength_results = run_pnfam_str_calculations()
    log("")
    log("Tests: {}, Passed: {}, Errors {}".format(tests, tests-errors, errors))
    log("Tolerance was {}".format(STRENGTH_TOLERANCE))
    log("")

    
    # Rotational checks
    log("Running STR rotational analysis")
    log("")
    tests, errors = run_pnfam_str_rot_analysis(strength_results)
    log("")
    log("Tests: {}, Passed: {}, Errors {}".format(tests, tests-errors, errors))
    log("Tolerance was {}".format(STRENGTH_TOLERANCE))
    log("")
