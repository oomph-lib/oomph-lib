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

# Validation for adaptive Helmholtz point source
#-----------------------------------------------
cd Validation

echo "Running adaptive Helmholtz point source"
mkdir RESLT
../adaptive_helmholtz_point_source > OUTPUT_adapt
echo "done"
echo " " >> validation.log
echo "Adpative Helmholtz point source" >> validation.log
echo "-------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > adaptive_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/adaptive_results.dat.gz   \
    adaptive_results.dat 0.5 1.0e-14 >> validation.log
fi
mv RESLT RESLT_adapt


# Validation for Helmholtz point source
#--------------------------------------

echo "Running Helmholtz point source"
mkdir RESLT
../helmholtz_point_source > OUTPUT_non_adapt
echo "done"
echo " " >> validation.log
echo "Helmholtz point source" >> validation.log
echo "----------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log
cat RESLT/trace.dat > non_adaptive_results.dat

if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/non_adaptive_results.dat.gz   \
    non_adaptive_results.dat  0.5 1.0e-14 >> validation.log
fi
mv RESLT RESLT_non_adapt



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
