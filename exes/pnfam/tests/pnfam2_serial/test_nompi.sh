#!/bin/bash

# An regression test script that solves a small-scale test problem
# In a realistic calculation, a model space of 16-20 shells should be used, but here
# we only use 6; This script is only meant for testing an installation on a new computer


# --------------------------------------------------------------------
#                           User Inputs
# --------------------------------------------------------------------
hfbtho_path=../../../hfbtho_gitlab/src/hfbtho
pnfam_path=../..

# --------------------------------------------------------------------
#                               Names
# --------------------------------------------------------------------
hfbtho_exe=hfbtho_main
pnfam_exe=contour_main.x
hfbtho_binary_output=hfbtho_output.hel
test_output_hfbtho=test_output_hfbtho.out
test_output_pnfam=test_output_pnfam.out
test_input_str=pnfam_CONTOUR.dat

# --------------------------------------------------------------------
#                            Test functions
# --------------------------------------------------------------------
test_file_exists () {
   # Checks if the file exists.
   # Arguments:
   # $1: file name (incl. path)
   # $2: message to show if the file doesn't exist (troubleshooting advice)
   echo -n "Checking that" $1 "exists... "
   if [ -f $1 ]
   then
      echo "OK"
      let tests_passed+=1
   else
      echo "FAILED"
      echo $2
      let tests_failed+=1
   fi
}

test_run_program () {
   # Runs a command/program, and checks if it completed without an error.
   # Arguments:
   # $1: description of the test
   # $2: the command to run
   # $3: message to show if the file doesn't exist (troubleshooting advice)
   echo -n $1"... "
   eval $2
   if [ $? -eq 0 ]
   then
      echo "OK"
      let tests_passed+=1
   else
      echo "FAILED"
      echo $3
      let tests_failed+=1
   fi
}

test_approximately_equal () {
   # Checks that the argument $1 is within 0.1% of $2; if not, $3 is printed as
   # a troubleshooting advice
   x=$(awk '{aux = ($1/($2) - 1); if (aux < 0) aux=(-aux); print (aux > 0.001)}' <<< "$1 $2")
   if [ $x -eq 0 ]
      then
      echo "OK"
      let tests_passed+=1
   else
      echo "FAILED"
      echo $3
      let tests_failed+=1
   fi
}

get_nth_word () {
   arr=($(echo $1))
   echo "${arr[$2]}"
}

get_nth_last_word () {
   arr=($(echo $1))
   let e=${#arr[@]}-$2
   echo "${arr[$e]}"
}

# --------------------------------------------------------------------
#                            Begin Testing
# --------------------------------------------------------------------
echo "----------------------------------------------------------------------"
echo "                  pnFAM (Serial) v2.0 regression test"
echo ""
echo " * This should take no more than a couple of minutes to complete."
echo " * Test outputs are piped to files 'test_output_*.out'."
echo "----------------------------------------------------------------------"
echo ""

let tests_passed=0
let tests_failed=0

# --------------------------------------------------------------------
#                          Executables Found?
# --------------------------------------------------------------------
test_file_exists "$hfbtho_path/$hfbtho_exe" "Check that $hfbtho_exe exists and that hfbtho_path is correctly set in this script."
test_file_exists "$pnfam_path/$pnfam_exe" "Check that $pnfam_exe exists and that pnfam_path is correctly set in this script."

if [ $tests_passed -ne 2 ]
then
   echo " *** Not all executables found. Aborting testing. ***"
   exit 1
fi

# --------------------------------------------------------------------
#                            HFBTHO Run
# --------------------------------------------------------------------
# Run the HFBTHO (the input parameters are in the hfbtho_NAMELIST.dat of the workind dir)
test_run_program "Running HFBTHO (this should take less than a minute)" \
                 "$hfbtho_path/$hfbtho_exe > $test_output_hfbtho" \
                 "There was an error running HFBTHO."

# Was hfbtho binary created? If not, stop the tests, the HFBTHO was probably not patched.
test_file_exists $hfbtho_binary_output "The output file $hfbtho_binary_output was not created by HFBTHO."

# Compare the total energy in the HFBTHO output to expected
# It's very unlikely you would get even close to the correct result by blind luck if
# there was a problem
echo -n "Checking the HFBTHO solution total energy matches expected... "
res=$(get_nth_last_word "$(grep --text 'tEnergy: ehfb (qp)...' thoout.dat)" 1)
test_approximately_equal $res "-432.187650" "The HFBTHO total energy was not within 0.1% of the expected value."

# --------------------------------------------------------------------
#                         PNFAM Strength Run
# --------------------------------------------------------------------
# Run the pnFAM in STR mode for a few points (the input parameters from default.in)
test_run_program "Running pnFAM strength (this should take a couple of minutes at most)" \
   "$pnfam_path/$pnfam_exe $test_input_str > $test_output_pnfam" \
   "There was an error running $pnfam_exe."

# Compare the pnFAM output to expected
echo -n "Checking pnFAM strength function value 1/3... "
line=$(tail -n 1 GT-K1.out)
val=$(get_nth_word "$line" 2)
test_approximately_equal $val "4.0138005453632780E-01" \
   "The pnFAM strength was not within 0.1% of the expected value."

echo -n "Checking pnFAM strength function value 2/3... "
line=$(tail -n 2 GT-K1.out | head -n 1)
val=$(get_nth_word "$line" 2)
test_approximately_equal $val "1.2528722692977165E-01" \
   "The pnFAM strength was not within 0.1% of the expected value."

echo -n "Checking pnFAM strength function value 3/3... "
line=$(tail -n 3 GT-K1.out | head -n 1)
val=$(get_nth_word "$line" 2)
test_approximately_equal $val "1.0604161621265991E-01" \
   "The pnFAM strength was not within 0.1% of the expected value."

# --------------------------------------------------------------------
#                            Test Results
# --------------------------------------------------------------------
let tests_total=tests_passed+tests_failed
echo ""
echo "Tests passed:" $tests_passed
echo "Tests failed:" $tests_failed

if [ $tests_failed -eq 0 ]
then
   echo ""
   echo " *** The serial pnFAM appears to be correctly compiled. ***"
   echo ""
fi
