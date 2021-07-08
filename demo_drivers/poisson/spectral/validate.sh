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

# Validation for demo poisson
#----------------------------
cd Validation

echo "Running 1D demo spectral poisson validation "
mkdir RESLT_1d
cd RESLT_1d
../../one_d_spectral > ../OUTPUT_one_d_spectral
cd ..
echo "done"
echo " " >> validation.log
echo "1D demo spectral poisson validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_1d/soln0.dat RESLT_1d/soln1.dat > 1d_spectral.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/1d_spectral.dat.gz   \
    1d_spectral.dat  >> validation.log
fi

echo "Running 2D demo spectral poisson validation "
mkdir RESLT_2d
../two_d_spectral_adapt > ./OUTPUT_two_d_spectral_adapt
echo "done"
echo " " >> validation.log
echo "2D demo spectral poisson validation" >> validation.log
echo "--------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_2d/soln0.dat > 2d_spectral.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/2d_spectral.dat.gz   \
    2d_spectral.dat  >> validation.log
fi

echo "Running eighth sphere spectral poisson validation "
mkdir RESLT_3d
../three_d_spectral_adapt lala > OUTPUT_three_d_spectral_adapt
echo "done"
echo " " >> validation.log
echo "Eighth sphere spectral poisson validation" >> validation.log
echo "--------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_3d/soln0.dat RESLT_3d/soln1.dat > 3d_spectral.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/3d_spectral.dat.gz   \
    3d_spectral.dat  0.1 1.0e-8 >> validation.log
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
