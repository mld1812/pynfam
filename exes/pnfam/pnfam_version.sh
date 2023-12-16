#------------------------------------------------------------------------------
# This bash script outputs the contents of a file that contains the fortran
# variable definitions for the pnfam version and git commit number.
#
# M.T. Mustonen and T. Shafer, UNC Chapel Hill, 2013-15
# E.M. Ney, UNC Chapel Hill, 2018-
#------------------------------------------------------------------------------

sha1=`git rev-parse --short HEAD`

date=`git log -1 --format="%ci" | cut -d" " -f1`

diff=$([ "`git diff --shortstat`" != "" ] && echo " + updates" || echo "")

cat <<HERE
! Version
character(len=*), parameter :: commit = "${sha1} [$date]$diff"
character(len=*), parameter :: version ='2.00'
HERE
