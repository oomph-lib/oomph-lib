#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=2

# Setup validation directory
#---------------------------
rm -rf Validation
mkdir Validation

#######################################################################

# Validation adaptive solid traction
#-----------------------------------

cd Validation

mkdir RESLT

echo "Running 3d adaptive solid traction "
../three_d_traction bla > OUTPUT_adaptive

echo "done"
echo " " >> validation.log
echo "3D Adaptive solid traction" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/traction1.dat RESLT/soln1.dat > result_adaptive.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_adaptive.dat.gz \
    result_adaptive.dat  0.1 1.0e-11 >> validation.log
fi

mv RESLT RESLT_adaptive

mkdir RESLT


echo "Running 3d nonadaptive solid traction "
../three_d_traction_non_adapt bla > OUTPUT_nonadaptive

echo "done"
echo " " >> validation.log
echo "3D Non-adaptive solid traction" >> validation.log
echo "------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/traction1.dat RESLT/soln1.dat > result_nonadaptive.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/result_nonadaptive.dat.gz \
    result_nonadaptive.dat  0.1  1.0e-11 >> validation.log
fi

mv RESLT RESLT_nonadaptive


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
