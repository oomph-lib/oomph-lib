#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)


#Set the number of tests to be checked
NUM_TESTS=2


# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

# Validation for demo linear_wave
#-------------------------
cd Validation

echo "Running 2D demo linear_wave validation "
mkdir RESLT_smooth
mkdir RESLT_impulsive

../two_d_linear_wave blabla > OUTPUT
echo "done"
echo " " >> validation.log
echo "2D demo linear_wave validation" >> validation.log
echo "-----------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT_smooth/soln0.dat RESLT_smooth/soln1.dat RESLT_smooth/soln2.dat \
RESLT_impulsive/soln0.dat RESLT_impulsive/soln1.dat RESLT_impulsive/soln2.dat \
 > results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz   \
    results.dat  5.1 1.0e-13 >> validation.log
fi

mv RESLT_smooth RESLT_smooth_two_d_linear_wave
mv RESLT_impulsive RESLT_impulsive_two_d_linear_wave

# Do flux validation
echo "Running 2D demo linear_wave flux validation "
mkdir RESLT
../two_d_linear_wave_flux la la la > OUTPUT
echo "done"
echo " " >> validation.log
echo "2D demo linear_wave flux validation" >> validation.log
echo "----------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/soln0.dat RESLT/soln1.dat RESLT/soln2.dat > results_flux.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results_flux.dat.gz   \
    results_flux.dat 0.1 1.0e-13  >> validation.log
fi

mv RESLT RESLT_two_d_linear_wave_flux

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
