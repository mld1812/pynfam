#!/usr/bin/env python
# Take a matrix QRPA output file and fold with a (default) 1-MeV Lorentzian.

import argparse
import numpy as np
import sys


parser = argparse.ArgumentParser(description="Fold matrix QRPA output with a Lorentzian")
parser.add_argument('--dx',                        type=float, default=0.1,  help='strength function sampling interval')
parser.add_argument('--min',                       type=float, default=0.0,  help='minimum EQRPA (default 0.0)')
parser.add_argument('--max',                       type=float, default=None, help='maximum EQRPA')
parser.add_argument('--square', '-s', action='store_true',                   help='square strength/matrix element entries')
parser.add_argument('--width',  '-w', metavar='W', type=float, default=1.0,  help='Lorentzian distribution width')
parser.add_argument('file',           metavar='FILE', type=str)
CLI = parser.parse_args()


data = np.loadtxt(CLI.file)
if len(data.shape) < 2:
   data = np.array([data])
if CLI.square:
   data[:,1:] = data[:,1:]**2


L_BOUND = CLI.min
if CLI.max is not None:
   R_BOUND = CLI.max
else:
   R_BOUND = np.ceil(data[:,0].max() + 3*CLI.width)


nintv  = np.ceil((R_BOUND-L_BOUND)/CLI.dx)
points = np.linspace(L_BOUND, R_BOUND, int(nintv+1), endpoint=True)

FMT_EQRPA = '{:>15.9f}'
FMT_STR   = '{:>18.9e}'


# Main loop
ncols = data.shape[1]
nrows = data.shape[0]

for pt in points:
   sys.stdout.write(FMT_EQRPA.format(pt))
   lor = [(CLI.width/2.0/np.pi) / ((pt-en)**2 + CLI.width**2/4.0) for en in data[:,0]]
   for ic in xrange(1, ncols):
      sys.stdout.write(FMT_STR.format(sum(data[:,ic]*lor)))
   sys.stdout.write('\n')
