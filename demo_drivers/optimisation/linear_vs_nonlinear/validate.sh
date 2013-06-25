#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=3


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation



# Validation for linear vs nonlinear optimisation
#-------------------------------------------------
cd Validation


echo "Running linear vs nonlinear optimisation demo "
mkdir RESLT
echo "S" >> tmp_input.dat
../two_d_poisson 3 < tmp_input.dat > OUTPUT_two_d_poisson
rm tmp_input.dat
echo "done"
echo " " >> validation.log
echo "Linear vs nonlinear optimisation demo" >> validation.log
echo "-------------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.cpp_style.dat  RESLT/soln1.cpp_style.dat  \
    > results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
echo "Compare results: " >> validation.log
../../../../bin/fpdiff.py ../validata/results.dat.gz   \
    results.dat  >> validation.log
echo "Compare C and C++ output: " >> validation.log
../../../../bin/fpdiff.py RESLT/soln0.cpp_style.dat RESLT/soln0.c_style.dat\
    >> validation.log
../../../../bin/fpdiff.py RESLT/soln1.cpp_style.dat RESLT/soln1.c_style.dat\
    >> validation.log
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
