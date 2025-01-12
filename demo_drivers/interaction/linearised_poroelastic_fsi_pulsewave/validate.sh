#! /bin/sh

# Get the OOPMH-LIB root directory from a makefile
OOMPH_ROOT_DIR=$(make -s --no-print-directory print-top_builddir)

#Set the number of tests to be checked
NUM_TESTS=1

# Setup validation directory
#---------------------------
touch Validation
rm -r -f Validation
mkdir Validation

cd Validation

# Validation for poroelastic linearised FSI pulswave
#---------------------------------------------------

echo "Running porelastic linearised FSI pulswave"
mkdir RESLT

# Run in validation mode
../linearised_poroelastic_fsi_pulsewave --validation > OUTPUT
echo "done"
echo " " >> validation.log
echo "Linearised poroelastic FSI pulswave" >> validation.log
echo "-----------------------------------" >> validation.log
echo " " >> validation.log
echo "Validation directory: " >> validation.log
echo " " >> validation.log
echo "  " `pwd` >> validation.log
echo " " >> validation.log


cat RESLT/regular_poro1.dat | cut -d ' ' -f 1-4 > poro1.dat
cat RESLT/regular_poro3.dat | cut -d ' ' -f 1-4 > poro3.dat

#cat RESLT/regular_fluid1.dat RESLT/regular_fluid3.dat RESLT/regular_poro1.dat RESLT/regular_poro3.dat > results.dat
cat RESLT/regular_fluid1.dat RESLT/regular_fluid3.dat poro1.dat poro3.dat > results.dat


if test "$1" = "no_fpdiff"; then
  echo "dummy [OK] -- Can't run fpdiff.py because we don't have python or validata" >> validation.log
else
../../../../bin/fpdiff.py ../validata/results.dat.gz   \
         results.dat  0.5 3.0e-4 >> validation.log
fi

# Append log to main validation log
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
