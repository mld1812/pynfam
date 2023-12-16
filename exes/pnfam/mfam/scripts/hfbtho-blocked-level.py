#!/usr/bin/env python

# From a search of blocking candidates, report the orbital which corresponds
# to the most-bound solution by scraping an HFBTHO output file.

import argparse
import re
import subprocess as sp

parser = argparse.ArgumentParser()
parser.add_argument('file', type=str, nargs='+')
parser.add_argument('--verbose', '-v', action='store_true')
CLI = parser.parse_args()

for ii, file in enumerate(CLI.file):
   if ii > 0: print
   print 'FILE: %s' % file

   # Find solutions
   cmd = sp.check_output('grep -n "tEnergy" "%s"' % file, shell=True)
   solutions = []

   for line in cmd.strip().split('\n'):
      pieces   = line.strip().split()
      line_num = int(pieces[0].strip(':'))
      energy   = float(pieces[-1])
      solutions.append((line_num, energy))

   if CLI.verbose:
      print 'SOLUTIONS FOUND:'
      print '\n'.join(['{:>14.6f} MeV'.format(s[1]) for s in solutions])

   # Find the orbital for the minmial solution
   minimum = sorted(solutions[1:], key=lambda s:s[1])[0]
   prev_id = solutions.index(minimum)-1

   if CLI.verbose:
      print
      print 'MINIMUM: {:>.6f} MeV'.format(minimum[1])

   # Find the blocked orbital in a bad way... overwrite the value every time
   blocked = None
   with open(file) as f:
      for i, line in enumerate(f):
         if i <= prev_id or i > minimum[0]:
            continue

         search = re.search(r'blocking: block=([0-9\s]+) state=([0-9\s]+).*([0-9]+)([+-])\[', line, re.I)
         if search is None:
            continue
         else:
            blocked = '{:d}{:s} [block={:d}, state={:d}]'.format(
               int(search.group(3)), search.group(4), int(search.group(1)),
               int(search.group(2)))

   print 'ORBITAL: {:s}'.format(blocked)
