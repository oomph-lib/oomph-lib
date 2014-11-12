#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

set -o errexit
no_fpdiff=false # `test "$1" = "no_fpdiff"`
set -o nounset

#Set the number of tests to be checked
NUM_TESTS=3


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo unsteady heat
#----------------------------------
cd Validation


# # Midpoint safe elements with bdf2
# # ============================================================
# echo "Running 2D unsteady heat midpoint elements with bdf2 time stepper validation"
# mkdir RESLT
# ../two_d_unsteady_heat_midpoint "bdf2" > OUTPUT_bdf2
# echo "done"
# echo >> validation.log
# echo "2D unsteady heat midpoint elements with bdf2 time stepper validation" >> validation.log
# echo "------------------------------------------------" >> validation.log
# echo >> validation.log
# echo "Validation directory: " >> validation.log
# echo >> validation.log
# echo "  " `pwd` >> validation.log
# echo >> validation.log
# cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat > result_bdf2.dat
# mv RESLT RESLT_bdf2

# if $no_fpdiff; then
#   echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
# else
# ../../../../bin/fpdiff.py ../validata/result.dat.gz \
#     result_bdf2.dat  >> validation.log
# fi


# # With midpoint
# # ============================================================
# echo "Running unsteady heat midpoint elements with midpoint time stepper validation"
# mkdir RESLT
# ../two_d_unsteady_heat_midpoint "midpoint" > OUTPUT_midpoint
# echo "done"
# echo >> validation.log
# echo "2D unsteady heat midpoint elements with midpoint time stepper validation " >> validation.log
# echo "------------------------------------------------" >> validation.log
# echo >> validation.log
# echo "Validation directory: " >> validation.log
# echo >> validation.log
# echo "  " `pwd` >> validation.log
# echo >> validation.log
# cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat > result_midpoint.dat
# mv RESLT RESLT_midpoint

# if $no_fpdiff; then
#   echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
# else
# ../../../../bin/fpdiff.py ../validata/result.dat.gz \
#     result_midpoint.dat  >> validation.log
# fi

# # With adaptive midpoint
# # ============================================================
# echo "Running unsteady heat midpoint elements with adaptive midpoint time stepper validation"
# mkdir RESLT
# ../two_d_unsteady_heat_midpoint "adaptive-midpoint" > OUTPUT_adaptive_midpoint
# echo "done"
# echo >> validation.log
# echo "2D unsteady heat midpoint elements with adaptive midpoint time stepper validation " >> validation.log
# echo "------------------------------------------------" >> validation.log
# echo >> validation.log
# echo "Validation directory: " >> validation.log
# echo >> validation.log
# echo "  " `pwd` >> validation.log
# echo >> validation.log
# cat RESLT/soln*.dat > result_adaptive_midpoint.dat
# mv RESLT RESLT_adaptive_midpoint

# if $no_fpdiff; then
#   echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
# else
# ../../../../bin/fpdiff.py ../validata/result_adaptive_midpoint.dat.gz \
#     result_adaptive_midpoint.dat  >> validation.log
# fi


# With midpoint bdf1 implementation
# ============================================================
echo "Running unsteady heat midpoint elements with midpoint time stepper validation by bdf1 implementation validation"
mkdir RESLT
../two_d_unsteady_heat_midpoint "midpoint-bdf1" > OUTPUT_midpoint_bdf1
echo "done"
echo >> validation.log
echo "2D unsteady heat midpoint elements with midpoint time stepper by bdf1 implementation validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo >> validation.log
echo "Validation directory: " >> validation.log
echo >> validation.log
echo "  " `pwd` >> validation.log
echo >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat > result_midpoint_bdf1.dat
mv RESLT RESLT_midpoint_bdf1

if $no_fpdiff; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result.dat.gz \
    result_midpoint_bdf1.dat  >> validation.log
fi


# With adaptive midpoint bdf1 implementation
# ============================================================
echo "Running unsteady heat midpoint elements with adaptive midpoint by bdf1 implementation time stepper validation"
mkdir RESLT
../two_d_unsteady_heat_midpoint "adaptive-midpoint-bdf1" > OUTPUT_adaptive_midpoint_bdf1
echo "done"
echo >> validation.log
echo "2D unsteady heat midpoint elements with adaptive midpoint by bdf1 implementation time stepper validation " >> validation.log
echo "------------------------------------------------" >> validation.log
echo >> validation.log
echo "Validation directory: " >> validation.log
echo >> validation.log
echo "  " `pwd` >> validation.log
echo >> validation.log
cat RESLT/soln*.dat > result_adaptive_midpoint_bdf1.dat
mv RESLT RESLT_adaptive_midpoint_bdf1

if $no_fpdiff; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_adaptive_midpoint.dat.gz \
    result_adaptive_midpoint_bdf1.dat  >> validation.log
fi

# Append output to global validation log file
#--------------------------------------------
cat validation.log >> ../../../../validation.log


cd ..



#######################################################################


#Check that we get the correct number of OKs
# validate_ok_count will exit with status
# 0 if all tests has passed.
# 1 if some tests failed.
# 2 if there are more 'OK' than expected.
. $OOMPH_ROOT_DIR/bin/validate_ok_count

# Never get here
exit 10
